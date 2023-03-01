// todo:
// a széléről érkezzenek a neutrofilok
// AB térbeli, széléről
// neutrofil time step -- improve

// az antibiotikum túl gyors, túl erős -- realisztikusabbá tenni, sztochasztikusan marad pár egy darabig?
// efficacy-k szétszedése. pl, egyelőre egy azonos immune efficacy van a decay-ben is és a source-ban is
// immunválasz lassú lecsengése

// mit szeretnénk kihozni:
// 1) ha 2-3 napon belül kap antibiotikumot, és mást nem, akkor nem károsodik a porc jelentősen
// 2) ha az 5-6. napon jut csak antibiotikumhoz, és máshoz nem, akkor jelentősen károsodhat a porc
// 3) ha az 5.-6. napon jut csak antibiotikumhoz, viszont ekkor kortikoszteroidot is kap, akkor nagyobb eséllyel marad meg a porcborítás

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
import java.util.Random;

public class CorticosteroidsAntibiotics{

    public static void main(String[] args) {

        int y = 100, x = 100, visScale = 3;
        boolean isAntibiotics = true;
        boolean isCorticosteroids = true;

        GridWindow win = new GridWindow("Cellular state space, bacterial concentration, immune reaction.", x*4, y, visScale,true);
        OpenGL2DWindow neutrophilWindow = new OpenGL2DWindow("Neutrophils", 500, 500, x,y);

        NewExperiment experiment = new NewExperiment(x, y, visScale, new Rand(1), isAntibiotics, isCorticosteroids, 24*60, 0.5, 4);
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
    public int visScale = 3;
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
    public double stapyloReproductionRate = 5.5 * Math.pow(10,-3);
    public double staphyloDiffCoeff; // D_V [sigma^2 / min]
    public double antibioticsDiffCoeff;

    Drug antibiotics;
    Drug corticosteroid;

    public double immuneResponseDecay = 0.00005;
    public boolean isAntibiotics = true;
    public boolean isCorticosteroid = true;

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
        this.antibiotics = new Drug();
        this.corticosteroid = new Drug();
        this.numberOfTicksDrug = 5 * 24 * 60; // we administer antibiotics for 5 days, i.e. 5*24*60 minutes

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

        bacterialCon.Add(initialPlace, 20.0);
        bacterialCon.Add(initialPlace-1, 20.0);
        bacterialCon.Update();

    }

    public double RunExperiment(GridWindow win, OpenGL2DWindow neutrophilWindow){

        double[] cellCounts = CountCells();

        for (int tick = 0; tick < this.numberOfTicks; tick ++){

            if (neutrophilWindow.IsClosed()){
                break;
            }
            neutrophilWindow.TickPause(20);

            TimeStep(tick);
            DrawModel(win);
            neutrophilLayer.DrawModel(neutrophilWindow);

            if( tick > 0 && ( (tick % 60) == 0 ))
                win.ToPNG(outputDir + "hour" + tick/(60) + ".jpg");

            double totalBacterialCon = TotalBacterialCon();
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
            if (cell.CellType == 0){
                healthyCells += 1;
            }
            else if (cell.CellType == 1 ){
                damagedCells += 1;
            }
            else if (cell.CellType == 2){
                deadCells += 1;
            }
        }

        cellCount[0] = healthyCells;
        cellCount[1] = damagedCells;
        cellCount[2] = deadCells;

        return cellCount;
    }

    void TimeStepNeutrophils(){

        // decay of the immune response level
        double decayFactor = immuneResponseDecay + 0.01 * corticosteroid.Efficacy(corticosteroidCon);
        for(Neutrophil cell: neutrophilLayer){
            cell.NeutrophilDecay(decayFactor);
        }

        // neutrophils arrive when bacterial concentration is high
        double corticosteroidSuppression = (1 - corticosteroid.Efficacy(corticosteroidCon));
        double bacterialFactor = 1 / (1 + 1/(Math.pow(TotalBacterialCon(),2)));
        double arrivalFactor = bacterialFactor * corticosteroidSuppression;
        BacteriaCallsNeutrophils(arrivalFactor);

        // neutrophils arrive when neutrophils call them by signaling
        for(Neutrophil cell: neutrophilLayer){
            Rand rng = new Rand();
            if(rng.Double() < neutrophilLayer.signalingIntensity) {
                cell.NeutrophilsCallNeutrophils(neutrophilLayer.signalSuccess * arrivalFactor);
            }
        }

        // neutrohils move inside the joint
        for(Neutrophil cell: neutrophilLayer){
            cell.NeutrophilMove();
        }
    }

    public void BacteriaCallsNeutrophils(double arrivalFactor){
        for (int i = 0; i < arrivalFactor; i++) {
            int[] arrivalPlace = Neutrophil.Infiltration();
            neutrophilLayer.NewAgentPT(arrivalPlace[0], arrivalPlace[1]).Init();
        }
    }
    void TimeStepStaphylo(){
        // if bacteria dies out
        if (this.TotalBacterialCon() < Math.pow(10, -6)){
            for (CartilageCell cell : cartilageLayer){
                bacterialCon.Set(cell.Isq(), 0);
            }
        }
        // bacterial reproduction
        for (CartilageCell cell : cartilageLayer){
            bacterialCon.Add(cell.Isq(), stapyloReproductionRate * bacterialCon.Get(cell.Isq()));
        }

        // diffusion of the bacteria
        bacterialCon.DiffusionADI(staphyloDiffCoeff);
        bacterialCon.Update();

        // removal of bacteria by antibiotics and immune cells
        for (CartilageCell cell : cartilageLayer){
            double immuneBacterialRemovalEff =  0.5 * 1 / (1 + 1/(Math.pow(neutrophilLayer.PopAt(cell.Isq()), 2)));
            double drugBacterialRemovalEff =  0.05 * antibiotics.Efficacy(antibioticsLayer.Get(cell.Isq()));
            System.out.println(drugBacterialRemovalEff);
            bacterialCon.Add(cell.Isq(), -drugBacterialRemovalEff * bacterialCon.Get(cell.Isq()));
            bacterialCon.Add(cell.Isq(), -immuneBacterialRemovalEff * bacterialCon.Get(cell.Isq()));
        }
        bacterialCon.Update();

    }

    void TimeStepDrug(int tick){

        // decay of the drugs
        this.antibioticsCon -= this.antibiotics.antibioticsDecay * this.antibioticsCon;
        this.corticosteroidCon -= this.corticosteroid.corticosteroidDecay * this.corticosteroidCon;

        // decay of the antibiotics in the first compartment
        // appearance of the antibiotics at the cartilage cells
        double transferQuantity = this.antibiotics.antibioticsDecayCompartment1 * this.antibioticsConCompartment1;
        this.antibioticsCon += transferQuantity;
        this.antibioticsConCompartment1 -= transferQuantity;

        // decay of the corticosteroids in the first compartment
        // appearance of the corticosteroids at the cartilage cells
        transferQuantity = this.corticosteroid.corticosteroidDecayCompartment1 * this.corticosteroidConCompartment1;
        this.corticosteroidCon += transferQuantity;
        this.corticosteroidConCompartment1 -= transferQuantity;

        // drug appearance in the first compartment
        this.antibioticsConCompartment1 += AntibioticsSourceCompartment1(tick);
        this.corticosteroidConCompartment1 += CorticosteroidSourceCompartment1(tick);

        // antibiotics infiltration through borders
        for (int i = 0; i < cartilageLayer.border.length; i++){
            antibioticsLayer.Set(cartilageLayer.border[i], this.antibioticsCon);
        }
        // local decay of the antibiotics
        for (CartilageCell cell : cartilageLayer){
            antibioticsLayer.Add(cell.Isq(), -this.antibiotics.antibioticsDecay * antibioticsLayer.Get(cell.Isq()));
        }
        antibioticsLayer.Update();

        // diffusion of the antibiotics
        antibioticsLayer.DiffusionADI(antibioticsDiffCoeff);
        antibioticsLayer.Update();

    }

    void TimeStepCartilageCells(){

        // immune cells damage healthy cells
        for (CartilageCell cell : cartilageLayer){
            cell.CellDamage(neutrophilLayer);
        }
        // damaged cartilage cells die
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
        for (int i = 0; i < cartilageLayer.length; i++){
            cellularBacterialCon[i] = bacterialCon.Get(i);
        }
        for (double bacterialConInCell : cellularBacterialCon ){
            totalBacterialCon = totalBacterialCon + bacterialConInCell;
        }
        return totalBacterialCon;
    }

    double AntibioticsSourceCompartment1(int tick){

        if ((tick > numberOfTicksDelay) && (isAntibiotics == true) && (((tick - numberOfTicksDelay) % (12 * 60)) == 1)) {
            return this.antibiotics.antibioticsSourceCompartment1;
        } else {
            return 0.0;
        }
    }

    double CorticosteroidSourceCompartment1(int tick){

        if ((tick > numberOfTicksDelay) && (isCorticosteroid == true) && (((tick - numberOfTicksDelay) % (12 * 60)) == 1)) {
            return this.corticosteroid.corticosteroidSourceCompartment1;
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
        String projPath = PWD() + "/output/StaphyloExperiments";
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
            drugInfo = this.antibiotics.antibioticsSourceCompartment1;
        }

        String outputDir = projPath + "/" + date_time + drugInfo + "__diff" + this.staphyloDiffCoeff;
        outputDir +=  "__delaytime" + this.numberOfTicksDelay + "/";

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
                if (drawMe.CellType == 0){		// Healthy cells
                    vis.SetPix(i, RGB256(119, 198, 110));
                }
                else if (drawMe.CellType == 1){   // Damaged cells
                    vis.SetPix(i, RGB256(124, 65, 120));
                }
                else if (drawMe.CellType == 2){
                    vis.SetPix(i, RGB(0, 0, 0));
                }
                else if (drawMe.CellType == 3){
                    vis.SetPix(i, RGB(255, 255, 255));
                }
            }

            vis.SetPix(cartilageLayer.ItoX(i) + x, cartilageLayer.ItoY(i), HeatMapGRB(100*bacterialCon.Get(i)));

            vis.SetPix(cartilageLayer.ItoX(i) + 2 * x, cartilageLayer.ItoY(i), HeatMapBGR(0.1*neutrophilLayer.PopAt(i)));

            vis.SetPix(cartilageLayer.ItoX(i) + 3 * x, cartilageLayer.ItoY(i), HeatMapGRB(0.1*antibioticsLayer.Get(i)));
        }
    }
}
class Drug{

    // general properties
    double EC50 = 20; // the drug quantity which results in 50% efficacy, in nM = nanoMolars, https://www.fda.gov/media/155050/download

    double corticosteroidDecay = 0.010;
    double corticosteroidSourceCompartment1 = 200;
    double corticosteroidDecayCompartment1 = 0.010;

    double antibioticsDecay = 0.025;
    double antibioticsSourceCompartment1 = 200;
    double antibioticsDecayCompartment1 = 0.010;


    double Efficacy(double drugNow){

        // 1 / (1 + 1/(Math.pow(antibioticsCon / 100, 2)));
        double Efficacy = 1 / ( 1 + (EC50 / StochasticDrug(drugNow)));
        return Efficacy;

    }

    double StochasticDrug(double drug){

        double stdDevOfGaussian = drug / 100;

        Rand random = new Rand();
        double stochasticDrug = random.Gaussian(drug, stdDevOfGaussian);
        stochasticDrug = stochasticDrug > 0.0 ? stochasticDrug : 0.0;

        return stochasticDrug;

    }

}
class CartilageLayer extends AgentGrid2D<CartilageCell> {
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
    int CellType;
    Rand rn = new Rand();
    double damageRate = Math.pow(10,-3);
    double deathProb = Math.pow(10,-4);

    public void CellInit(boolean isHealthy, boolean isDamaged,
                         boolean isDead){

        if(isHealthy == true){
            this.CellType = 0;
        }
        else if(isDamaged == true){
            this.CellType = 1;
        }
        else if(isDead == true) {
            this.CellType = 2;
        }
    }

    public void CellDamage(NeutrophilLayer neutrophilLayer){

        // nor antibiotics nor corticosteroids help the cells directly
        // difference w.r.t. certain antivirals, which prevent the cells from getting infected

        // chondral damage: depends on the immune response strength (a strong immune response kills bacteria, but damages the cartilage as well)

        double immuneConAtCell = neutrophilLayer.PopAt(Isq());
        double effectiveDamageProb = immuneConAtCell * damageRate;

        if (this.CellType == 0){ // healthy cell
            if (rn.Double() < effectiveDamageProb) {
                this.CellType = 1;
            }
        }
    }

    public void CellDeath(){

        if (this.CellType == 1) { // damaged
            if(rn.Double() < deathProb){
                this.CellType = 2;
            }
        }
    }
}

class Neutrophil extends AgentPT2D<NeutrophilLayer> {
    int color;
    Rand random = new Rand();

    public void Init() {this.color = RGB256(random.Int(255),random.Int(255),random.Int(255));}

    public void NeutrophilDecay(double decayProb){
        if(G.rng.Double()<decayProb){
            Dispose();
        }
    }

    public void NeutrophilsCallNeutrophils(double signalSuccess){
        if((G.PopAt(Isq())<1) && G.rng.Double()<signalSuccess){ //
            int[] arrivalPlace = Infiltration();
            G.NewAgentPT(arrivalPlace[0], arrivalPlace[1]).Init();
        }
    }

    public void NeutrophilMove(){
        G.rng.RandomPointInCircle(2,G.moveCoords);
        MoveSafePT(Xpt()+G.moveCoords[0], Ypt()+G.moveCoords[1]);
    }

    public static int[] Infiltration(){
        Rand random = new Rand();
        int randomX, randomY;
        int randomNumber = random.Int(4);
        if(randomNumber < 1){
            randomX = 1;
            randomY = random.Int(100);
        }
        else if(randomNumber < 2){
            randomX = random.Int(100);
            randomY = 1;
        }
        else if(randomNumber < 3){
            randomX = random.Int(100);
            randomY= 99;
        }
        else{
            randomX = 99;
            randomY = random.Int(100);
        }
        int[] arrivalPlace = {randomX, randomY};
        return arrivalPlace;
    }

}

class NeutrophilLayer extends AgentGrid2D<Neutrophil> {

    double signalingIntensity = 0.5;
    double signalSuccess = 0.7;

    Rand rng = new Rand();
    double[]moveCoords=new double[2];

    public NeutrophilLayer(int x, int y) {super(x, y, Neutrophil.class);}

    public void DrawModel(OpenGL2DWindow win){
        win.Clear(Util.BLACK);
        for (Neutrophil cell : this){
            win.Circle(cell.Xpt(), cell.Ypt(),0.5, cell.color);
        }
        win.Update();
    }

}