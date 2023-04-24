package CorticosteroidsAntibiotics;

import HAL.GridsAndAgents.AgentSQ2Dunstackable;
import HAL.GridsAndAgents.AgentGrid2D;
import HAL.GridsAndAgents.PDEGrid2D;
import HAL.Gui.GridWindow;
import HAL.Tools.FileIO;
import HAL.Rand;
import static HAL.Util.*;

import HAL.GridsAndAgents.AgentPT2D;
import HAL.Gui.OpenGL2DWindow;
import HAL.Util;

import static HAL.Util.RGB256;

import java.io.File;

public class CorticosteroidsAntibiotics{

    public static void main(String[] args) {

        int y = 100, x = 100, visScale = 3;
        boolean isAntibiotics = true;
        boolean isCorticosteroids = true;

        GridWindow win = new GridWindow("Cellular state space, bacterial concentration, immune reaction.", x*4, y, visScale,true);
        OpenGL2DWindow neutrophilWindow = new OpenGL2DWindow("Neutrophils", 500, 500, x,y);

        NewExperiment experiment = new NewExperiment(x, y, visScale, new Rand(1), isAntibiotics, isCorticosteroids, 150, 0.6, 20);
        experiment.numberOfTicks = experiment.numberOfTicksDelay + experiment.numberOfTicksDrug;

        experiment.Init();
        double remainingHealthyCells = experiment.RunExperiment(win, neutrophilWindow);

        System.out.println("Remaining healthy cells: " + remainingHealthyCells);

        win.Close();

    }

}

class NewExperiment{

    public int x = 100;
    public int y = 100;
    public int visScale;
    public int numberOfTicksDelay;
    public int numberOfTicksDrug;
    public int numberOfTicks;
    public PDEGrid2D bacterialCon;
    public PDEGrid2D antibioticsLayer;
    public NeutrophilLayer neutrophilLayer;
    public CartilageLayer cartilageLayer;
    //public PDEGrid2D immuneResponseLevel; // similar to T-cell concentrations, but more generic
    public double antibioticsCon = 0;
    public double antibioticsConCompartment1 = 0;
    public double corticosteroidCon = 0;
    public double corticosteroidConCompartment1 = 0;
    public Rand rn;
    public double[] cellularBacterialCon = new double[x*y];
    public double staphyloReproductionRate = Math.pow(10,-2);
    public double staphyloDiffCoeff; // D_V [sigma^2 / min]
    public double antibioticsDiffCoeff;
    public boolean staphyloExtinct = false;
    public double extinctDeclarationThreshold = 0.05;

    Drug antibiotics;
    Drug corticosteroid;

    public double immuneResponseDecay = 0.0005;

    public boolean isAntibiotics;
    public boolean isCorticosteroid;

    public FileIO outFile;
    public FileIO paramFile;
    public FileIO concentrationsFile;
    public String outputDir;

    public NewExperiment(int xDim, int yDim, int visScale, Rand rn, boolean isAntibiotics, boolean isCorticosteroid, int numberOfTicksDelay, double staphyloDiffCoeff, double antibioticsDiffCoeff){

        this.x = xDim;
        this.y = yDim;
        this.visScale = visScale;
        this.numberOfTicksDelay = numberOfTicksDelay;
        this.staphyloDiffCoeff = staphyloDiffCoeff;
        this.antibioticsDiffCoeff = antibioticsDiffCoeff;
        this.rn = rn;

        this.isAntibiotics = isAntibiotics;
        this.isCorticosteroid = isCorticosteroid;

        this.antibiotics = new Drug("antibiotics", 0.025, 100, 0.010,20);
        this.corticosteroid = new Drug("corticosteroid", 0.010, 50, 0.010,20);
        this.numberOfTicksDrug = 5000;

        bacterialCon = new PDEGrid2D(xDim, yDim);
        bacterialCon.Update();

        antibioticsLayer = new PDEGrid2D(xDim, yDim);
        antibioticsLayer.Update();

        neutrophilLayer = new NeutrophilLayer(xDim, yDim);
        cartilageLayer = new CartilageLayer(xDim, yDim);

        outputDir = this.OutputDirectory();
        outFile = new FileIO(outputDir.concat("/").concat("Out").concat(".csv"), "w");
        paramFile = new FileIO(outputDir.concat("/").concat("Param").concat(".csv"), "w");
        concentrationsFile = new FileIO(outputDir.concat("/").concat("concentrations").concat(".csv"), "w");

    }

    public void Init(){

        WriteHeader();
        cartilageLayer.Init();
        InfectionInit();
    }

    public void InfectionInit(){

        int initialPlace = this.cartilageLayer.length / 2;

        bacterialCon.Add(initialPlace, 25.0);
        bacterialCon.Add(initialPlace-1, 2.0);
        bacterialCon.Update();

    }

    public double RunExperiment(GridWindow win, OpenGL2DWindow neutrophilWindow){

        double[] cellCounts = CountCells();
        // System.out.println(cellCounts[0]+", " + cellCounts[1] + ", " + cellCounts[2]);

        for (int tick = 0; tick < this.numberOfTicks; tick ++){

            if (neutrophilWindow.IsClosed()){
                break;
            }
            neutrophilWindow.TickPause(20);

            TimeStep(tick);
            DrawModel(win);
            neutrophilLayer.DrawModel(neutrophilWindow);

            if( tick > 0 && ( (tick % (1*10)) == 1 )) {
                win.ToPNG(outputDir + tick + ".jpg");
                neutrophilWindow.ToPNG(outputDir + tick + "neutrophils.jpg");
            }

            double totalBacterialCon = TotalBacterialCon();
            System.out.println("TBC: " + totalBacterialCon);
            double totalImmuneResponseLevel = neutrophilLayer.Pop();
            cellCounts = CountCells();
            concentrationsFile.Write(totalBacterialCon + "," + totalImmuneResponseLevel + "," + antibioticsCon + "," + corticosteroidCon + "\n");
            outFile.Write(tick +"," + cellCounts[0] + "," + cellCounts[1]+
                    "," + cellCounts[2] + "," + totalBacterialCon + "," + totalImmuneResponseLevel + "," + antibioticsCon + "," + corticosteroidCon + "\n");

        }

        CloseFiles();
        neutrophilWindow.Close();

        return cellCounts[0];

    }

    double[] CountCells(){

        double healthyCells = 0, damagedCells = 0, deadCells = 0;
        double[] cellCount = new double[3];
        for (CartilageCell cell: cartilageLayer){
            if (cell.cellType == 0){
                healthyCells += 1;
            }
            else if (cell.cellType == 1 ){
                damagedCells += 1;
            }
            else if (cell.cellType == 2){
                deadCells += 1;
            }
        }

        cellCount[0] = healthyCells;
        cellCount[1] = damagedCells;
        cellCount[2] = deadCells;

        return cellCount;
    }

    void TimeStepNeutrophils(){

        // decay of the immuneResponseLevel
        double decayFactor = immuneResponseDecay + 0.05 * corticosteroid.DrugEfficacy(corticosteroidCon);
        for(Neutrophil cell: neutrophilLayer){
            cell.NeutrophilDecay(decayFactor);
        }

        // neutrophils arrive when bacterial concentration is high or when other neutrophils call them
        double corticosteroidSuppression = (1 - corticosteroid.DrugEfficacy(corticosteroidCon));
        double bacterialFactor = 10 / (1 + 1/(Math.pow(TotalBacterialCon(),2)));
        double arrivalFactor = bacterialFactor * corticosteroidSuppression;

        neutrophilLayer.BacteriaGeneratedNeutrophilArrival(arrivalFactor);

        // neutrophils "call" new ones by signaling, so this could be reasonable
        for(Neutrophil cell: neutrophilLayer){
            Rand rng = new Rand();
            double randomDouble = rng.Double();
            if(randomDouble < neutrophilLayer.signalingIntensity) {
               cell.NeutrophilsCallNeutrophils(300*neutrophilLayer.signalSuccess * corticosteroidSuppression * bacterialCon.Get(cell.Isq()));
           }
        }

        for(Neutrophil cell: neutrophilLayer){
            cell.NeutrophilMove();
        }

    }

    void TimeStepStaphylo(){

        if(this.staphyloExtinct == false) {

            // if total bacterial concentration is extremely low, we assume the infection
            // dies out
            if (this.TotalBacterialCon() < this.extinctDeclarationThreshold) {
                for (CartilageCell cell : cartilageLayer) {
                    bacterialCon.Set(cell.Isq(), 0);
                }
                this.staphyloExtinct = true;
            }

            // bacterial reproduction
            for (CartilageCell cell : cartilageLayer) {
                bacterialCon.Add(cell.Isq(), staphyloReproductionRate * bacterialCon.Get(cell.Isq()));
            }

            bacterialCon.DiffusionADI(staphyloDiffCoeff);
            bacterialCon.Update();

            // removal of bacteria
            for (CartilageCell cell : cartilageLayer) {
                // double removalEfficacy = 2/(1+Math.exp(100*drugNow));
                // double removalEfficacy = 100*Math.pow(drugNow, 2)/(1+100*Math.pow(drugNow,2));
                double drugBacterialRemovalEff = 0.1 * antibiotics.DrugEfficacy(antibioticsLayer.Get(cell.Isq()));
                //System.out.println(drugBacterialRemovalEff);
                double immuneBacterialRemovalEff = 0.05 * 1 / (1 + 1 / (Math.pow(neutrophilLayer.PopAt(cell.Isq()), 2)));
                bacterialCon.Add(cell.Isq(), -drugBacterialRemovalEff * bacterialCon.Get(cell.Isq()));
                bacterialCon.Add(cell.Isq(), -immuneBacterialRemovalEff * bacterialCon.Get(cell.Isq()));
            }
            bacterialCon.Update();
        }

    }

    void TimeStepDrug(int tick){

        // STEP 1: CALCULATE SERUM CONCENTRATION LEVELS

        // decay of the drug
        this.antibioticsCon -= this.antibiotics.drugDecay * this.antibioticsCon;
        this.corticosteroidCon -= this.corticosteroid.drugDecay * this.corticosteroidCon;

        // decay of the drug in the first compartment (stomach or place of injection)
        // and appearance of the drug in the blood everywhere, including veins near the joint
        double transferQuantity = this.antibiotics.drugDecayCompartment1 * this.antibioticsConCompartment1;
        this.antibioticsCon += transferQuantity;
        this.antibioticsConCompartment1 -= transferQuantity;

        transferQuantity = this.corticosteroid.drugDecayCompartment1 * this.corticosteroidConCompartment1;
        this.corticosteroidCon += transferQuantity;
        this.corticosteroidConCompartment1 -= transferQuantity;

        // drug appearance in Compartment1
        this.antibioticsConCompartment1 += AntibioticsSourceCompartment1(tick);
        this.corticosteroidConCompartment1 += CorticosteroidSourceCompartment1(tick);

        // STEP 2: CALCULATE CONCENTRATION LEVELS WITHIN THE JOINT

        // dynamics at the border
        for (int i = 0; i < cartilageLayer.border.length; i++){
            // whatever serum concentration levels we have outide,
            // we assume we have that exact concentration level at the border cells
            antibioticsLayer.Set(cartilageLayer.border[i], this.antibioticsCon);
        }

        //decay

        for (CartilageCell cell : cartilageLayer){
            antibioticsLayer.Add(cell.Isq(), -this.antibiotics.drugDecay * antibioticsLayer.Get(cell.Isq()));
        }
        antibioticsLayer.Update();

        // diffusion
        antibioticsLayer.DiffusionADI(antibioticsDiffCoeff);
        antibioticsLayer.Update();

    }

    void TimeStepCartilageCells(){

        for (CartilageCell cell : cartilageLayer){
            cell.CellDamage(neutrophilLayer);
        }
        for (CartilageCell cell : cartilageLayer) {
            cell.CellDeath();
        }

    }

    void TimeStep(int tick){

        TimeStepStaphylo();
        TimeStepCartilageCells();
        TimeStepNeutrophils(); // immune response
        TimeStepDrug(tick); // both antibiotics and steroid

    }

    double TotalBacterialCon() {

        double totalBacterialCon = 0;
        for (int i =0; i < cartilageLayer.length; i++){
            cellularBacterialCon[i] = bacterialCon.Get(i);
        }
        for (double bacterialConInCell : cellularBacterialCon){
            totalBacterialCon = totalBacterialCon + bacterialConInCell;
        }
        return totalBacterialCon;
    }

    double AntibioticsSourceCompartment1(int tick){

        if ((tick > numberOfTicksDelay) && (isAntibiotics == true) && (((tick - numberOfTicksDelay) % 700) == 1)) {
            return this.antibiotics.drugSourceCompartment1;
        } else {
            return 0.0;
        }
    }

    double CorticosteroidSourceCompartment1(int tick){

        if ((tick > numberOfTicksDelay) && (isCorticosteroid == true) && (((tick - numberOfTicksDelay) % 700) == 1)) {
            return this.corticosteroid.drugSourceCompartment1;
        } else {
            return 0.0;
        }
    }

    void WriteHeader(){

        outFile.Write("tick, healthy cells, damaged cells, dead cells, "
                + "total bacterial conc., total immune response, total AB, total CS \n");
    }

    void CloseFiles(){

        outFile.Close();
        paramFile.Close();
        concentrationsFile.Close();

    }

    String OutputDirectory(){

        java.util.Date now = new java.util.Date();
        java.text.SimpleDateFormat dateFormat = new java.text.SimpleDateFormat("yyyy-MM-dd_HH-mm-ss");
        String date_time = dateFormat.format(now);
        String projPath = PWD() + "/CorticoSteroidsAntibiotics/output/StaphyloExperiments";
        if (this.isAntibiotics == false){
            projPath += "/noDrug";
        } else if (this.isCorticosteroid == true){
            projPath += "/steroidBoostedAB";
        } else {
            projPath += "/AntibioticsOnly";
        }

        double drugInfo;
        if (this.isAntibiotics == false){
            drugInfo = 0.0;
        } else {
            drugInfo = this.antibiotics.drugSourceCompartment1;
        }

        String outputDir = projPath + "/";

        new File(outputDir).mkdirs();

        return outputDir;

    }

    void DrawModel(GridWindow vis){

        for (int i = 0; i < cartilageLayer.length; i++) {
            CartilageCell drawMe = cartilageLayer.GetAgent(i);

            if (drawMe == null){
                vis.SetPix(i, RGB256(255, 255, 255));
            }
            else{
                if (drawMe.cellType == 0){		// Healthy cells
                    vis.SetPix(i, RGB256(119, 198, 110));
                }
                else if (drawMe.cellType == 1){   // Damaged cells
                    vis.SetPix(i, RGB256(124, 65, 120));
                }
                else if (drawMe.cellType == 2){
                    vis.SetPix(i, RGB(0, 0, 0));
                }
                else if (drawMe.cellType == 3){
                    vis.SetPix(i, RGB(255, 255, 255));
                }
            }

            vis.SetPix(cartilageLayer.ItoX(i) + x, cartilageLayer.ItoY(i), HeatMapGRB(100*bacterialCon.Get(i)));

            vis.SetPix(cartilageLayer.ItoX(i) + 2 * x, cartilageLayer.ItoY(i), HeatMapBGR(0.1*neutrophilLayer.PopAt(i)));

            vis.SetPix(cartilageLayer.ItoX(i) + 3 * x, cartilageLayer.ItoY(i), HeatMapRBG(0.1*antibioticsLayer.Get(i)));
        }
    }

}

class CartilageLayer extends AgentGrid2D<CartilageCell>{

    public int[] border;

    public CartilageLayer(int xDim, int yDim){
        super(xDim, yDim, CartilageCell.class);
        border = new int[2*(xDim+yDim)];
    }

    public void Init(){

        int borderIndex = 0;
        for (int i = 0; i < this.length; i++){
            CartilageCell cartilageCell = NewAgentSQ(i);
            cartilageCell.CellInit(true,false, false);
            if (cartilageCell.Xsq() == 0 || cartilageCell.Xsq() == 99 || cartilageCell.Ysq() == 0 || cartilageCell.Ysq() == 99){
                border[borderIndex] = i;
                borderIndex += 1;
            }
        }

    }

}

class CartilageCell extends AgentSQ2Dunstackable<CartilageLayer>{
    int cellType;
    Rand rn = new Rand();
    double damageRate = Math.pow(10,-3);
    double deathProb = Math.pow(10,-4);

    public void CellInit(boolean isHealthy, boolean isDamaged,
                         boolean isDead){

        if(isHealthy == true){
            this.cellType = 0;
        }
        else if(isDamaged == true){
            this.cellType = 1;
        }
        else if(isDead == true) {
            this.cellType = 2;
        }
    }

    public void CellDamage(NeutrophilLayer neutrophilLayer){

        // nor antibiotics nor corticosteroids help the cells directly
        // difference w.r.t. certain antivirals, which prevent the cells from getting infected

        // chondral damage: depends on the immune response strength (a strong immune response kills bacteria, but damages the cartilage as well)

        double immuneConAtCell = neutrophilLayer.PopAt(Isq());
        double effectiveDamageProb = immuneConAtCell * damageRate;

        if (this.cellType == 0){ // healthy cell
            if (rn.Double() < effectiveDamageProb) { //rn.Double() returns a value from 0 to 1
                this.cellType = 1;
            }
        }
    }

    public void CellDeath(){

        if (this.cellType == 1) { // damaged
            if(rn.Double() < deathProb){
                this.cellType = 2;
            }
        }
    }
}

class Drug{

    // general properties
    String drugType;
    double EC50; // the drug quantity which results in 50% efficacy

    double drugDecay;
    double drugSourceCompartment1;
    double drugDecayCompartment1;

    public Drug(String drugType, double drugDecay, double drugSourceCompartment1, double drugDecayCompartment1, double EC50){

        this.drugType = drugType;
        this.drugDecay = drugDecay;
        this.drugSourceCompartment1 = drugSourceCompartment1;
        this.drugDecayCompartment1 = drugDecayCompartment1;
        this.EC50 = EC50;

    }

    public double DrugEfficacy(double drugNow){

        double drugEfficacy = 1 / ( 1 + (EC50 / StochasticDrug(drugNow)));
        return drugEfficacy;

    }

    double StochasticDrug(double drug){

        double stdDevOfGaussian = drug / 100;

        Rand random = new Rand();
        double stochasticDrug = random.Gaussian(drug, stdDevOfGaussian);
        stochasticDrug = stochasticDrug > 0.0 ? stochasticDrug : 0.0;

        return stochasticDrug;

    }

}

class Neutrophil extends AgentPT2D<NeutrophilLayer> {
    int color;
    Rand random = new Rand();

    public void Init() {this.color = RGB256(random.Int(255),random.Int(255),random.Int(255));}

    public void NeutrophilDecay(double decayProb){
        if(G.rn.Double()<decayProb){
            Dispose();
        }
    }

    public void NeutrophilsCallNeutrophils(double signalSuccess){
        if((G.PopAt(Isq())<5) && G.rn.Double()<signalSuccess){
            double[] location = G.NeutrophilInfiltrationLocation();
            G.NewAgentPT(location[0], location[1]).Init();
        }
    }

    public void NeutrophilMove(){
        G.rn.RandomPointInCircle(1.5,G.moveCoords);
        MoveSafePT(Xpt()+G.moveCoords[0], Ypt()+G.moveCoords[1]);
    }

}

class NeutrophilLayer extends AgentGrid2D<Neutrophil> {

    double signalingIntensity = 0.5;
    double signalSuccess = 0.02;

    Rand rn = new Rand();
    double[] moveCoords = new double[2];

    public NeutrophilLayer(int x, int y) {super(x, y, Neutrophil.class);}

    public void DrawModel(OpenGL2DWindow win){
        win.Clear(Util.BLACK);
        for (Neutrophil cell : this){
            win.Circle(cell.Xpt(), cell.Ypt(),0.5, cell.color);
        }
        win.Update();
    }

    public void BacteriaGeneratedNeutrophilArrival(double arrivalFactor){

        if(arrivalFactor > 1.0) {
            for (int i = 0; i < arrivalFactor; i++) {
                double[] location = this.NeutrophilInfiltrationLocation();
                this.NewAgentPT(location[0], location[1]).Init();
            }
        }

    }

    public double[] NeutrophilInfiltrationLocation(){

        double[] neutrophilRandLocation = new double[2];
        Rand random = new Rand();
        int edge = random.Int(4);

        if (edge == 0){
            neutrophilRandLocation[0] = random.Double(99);
            neutrophilRandLocation[1] = 99.0;
        } else if (edge == 1) {
            neutrophilRandLocation[0] = random.Double(99);
            neutrophilRandLocation[1] = 0;
        } else if (edge == 2) {
            neutrophilRandLocation[0] = 99;
            neutrophilRandLocation[1] = random.Double(99);
        } else {
            neutrophilRandLocation[0] = 0;
            neutrophilRandLocation[1] = random.Double(99);
        }

        return neutrophilRandLocation;
    }

}