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

    public static final int BIG_VALUE = (Integer.MAX_VALUE-1) / 2;

    public static void main(String[] args) {

        String singularOrSweep = "singular"; // "singular" or "sweep"
        String tabletsOrIntravenous = "tablets"; // drug taken orally or intravenously

        int y = 100, x = 100, visScale = 3;
        boolean isAntibiotics = true;
        boolean isCorticosteroids = true;

        GridWindow win = new GridWindow("Cellular state space, bacterial concentration, immune reaction.", x*3, y, visScale,true);
        OpenGL2DWindow neutrophilWindow = new OpenGL2DWindow("Neutrophils", 500, 500, x,y);

        if (singularOrSweep.equals("singular")){

            NewExperiment experiment = new NewExperiment(x, y, visScale, new Rand(1), isAntibiotics, isCorticosteroids, 1*6*60, 0.5, tabletsOrIntravenous, 110.0);
            experiment.numberOfTicks = experiment.numberOfTicksDelay + experiment.numberOfTicksDrug;

            experiment.Init();
            double remainingHealthyCells = experiment.RunExperiment(win, neutrophilWindow);

            if (isAntibiotics) {
                System.out.println(tabletsOrIntravenous.equals("tablets") ? "Tablets AB drug source [ng / ml]: " + experiment.antibiotics.drugSourceStomach : "Intravenous drug concentration [nM]: " + experiment.antibiotics.NgPerMlToNanomolars(experiment.antibiotics.intravenousAntibioticsCon));
            }
            if (isCorticosteroids) {
                System.out.println(tabletsOrIntravenous.equals("tablets") ? "Tablets CS drug source [ng / ml]: " + experiment.corticosteroid.drugSourceStomach : "Intravenous drug concentration [nM]: " + experiment.corticosteroid.NgPerMlToNanomolars(experiment.corticosteroid.intravenousCorticosteroidCon));
            }
            System.out.println("Remaining healthy cells: " + remainingHealthyCells);

        } else if (singularOrSweep.equals("sweep")) {

            String collectiveOutputDir = collectiveOutputDir(isAntibiotics, isCorticosteroids);
            FileIO collectiveOutFile = new FileIO(collectiveOutputDir.concat("/").concat("collectiveRemainingHealthyCells").concat(".csv"), "w");
            String collectiveResults;

            for (double staphyloDiffCoeffSweep = 0.00625; staphyloDiffCoeffSweep < 50; staphyloDiffCoeffSweep *= 2) {

                collectiveResults = "";
                collectiveResults += staphyloDiffCoeffSweep + ", ";

                for (double damageRateSweep = 0.0; damageRateSweep < 100.0; damageRateSweep += 5.0) {

                    NewExperiment experiment = new NewExperiment(x, y, visScale, new Rand(1), isAntibiotics, isCorticosteroids, BIG_VALUE, staphyloDiffCoeffSweep, tabletsOrIntravenous, damageRateSweep);
                    experiment.numberOfTicks = experiment.numberOfTicksDelay + experiment.numberOfTicksDrug;

                    experiment.Init();
                    double remainingHealthyCells = experiment.RunExperiment(win, neutrophilWindow);

                    collectiveResults += remainingHealthyCells + ", ";

                }
                collectiveOutFile.Write(collectiveResults + "\n");
                System.out.println(collectiveResults);
            }

            collectiveOutFile.Close();

        } else {
            System.out.println("The only options are singular and sweep.");
        }

        win.Close();

    }

    public static String collectiveOutputDir(boolean isAntibiotics, boolean isCorticosteroid){

        java.util.Date now = new java.util.Date();
        java.text.SimpleDateFormat dateFormat = new java.text.SimpleDateFormat("yyyy-MM-dd_HH-mm-ss");
        String date_time = dateFormat.format(now);
        String projPath = PWD() + "/output/StaphyloExperiments";
        if (isAntibiotics == false){
            projPath += "/noDrug";
        } else if (isCorticosteroid == true){
            projPath += "/steroidBoostedAB";
        } else {
            projPath += "/AntibioticsOnly";
        }
        String collectiveOutputDir = projPath + "/" + date_time + "__collective";
        new File(collectiveOutputDir).mkdirs();

        return collectiveOutputDir;

    }

}

class CorticosteroidDrug{

    // general properties
    double EC50 = 62; // in nM = nanoMolars, [nM] = 10^-9 [mol/L]; https://www.fda.gov/media/155050/download
    double molarMassDrug = 499.535;

    // intravenous properties
    double intravenousCorticosteroidCon = 5;

    // intravenous constructor
    public CorticosteroidDrug(double intravenousCorticosteroidCon, boolean isCorticosteroids){

        if (isCorticosteroids == true){
            this.intravenousCorticosteroidCon = intravenousCorticosteroidCon;
        } else {
            this.intravenousCorticosteroidCon = 0.0;
        }

    }

    // tablets case properties

    double drugDecay = 0.010;
    double drugSourceStomach = 400;
    double drugDecayStomach = 0.020;

    // tablets case constructor
    public CorticosteroidDrug(){

    }

    public double CorticosteroidEfficacy(double steroidNow){

        // 7000 * Math.pow(drugNow, 2)/(1+7000*Math.pow(drugNow,2));
        // drugNow is the drug concentration in nanograms / ml
        // drugNow needs to be converted to [nM]s, as IC50 is given in [nM]s
        double steroidNowInNanoMolars = NgPerMlToNanomolars(steroidNow);
        double steroidEfficacy = 1 / ( 1 + (EC50 / StochasticDrug(steroidNowInNanoMolars)));

        return steroidEfficacy;

    }

    double StochasticDrug(double drug){

        double stdDevOfGaussian = drug / 100;

        Rand random = new Rand();
        double stochasticDrug = random.Gaussian(drug, stdDevOfGaussian);
        stochasticDrug = stochasticDrug > 0.0 ? stochasticDrug : 0.0;

        return stochasticDrug;

    }

    double NgPerMlToNanomolars(double drugNow){
        return drugNow * Math.pow(10,3) / molarMassDrug;
    }

}

class AntibioticsDrug {

    // general properties
    double EC50 = 62; // in nM = nanoMolars, [nM] = 10^-9 [mol/L]; https://www.fda.gov/media/155050/download
    double molarMassDrug = 499.535;

    // intravenous properties
    double intravenousAntibioticsCon = 5;

    // intravenous constructor
    public AntibioticsDrug(double intravenousAntibioticsCon, boolean isAntibiotics){

        if (isAntibiotics == true){
            this.intravenousAntibioticsCon = intravenousAntibioticsCon;
        } else {
            this.intravenousAntibioticsCon = 0.0;
        }

    }

    // tablets decay properties

    double drugDecay = 0.005;
    double drugSourceStomach = 100;
    double drugDecayStomach = 0.010;

    // tablets case constructor
    public AntibioticsDrug(){

    }

    double AntibioticsEfficacy(double antibioticsNow){

        // 7000 * Math.pow(drugNow, 2)/(1+7000*Math.pow(drugNow,2));
        // drugNow is the drug concentration in nanograms / ml
        // drugNow needs to be converted to [nM]s, as IC50 is given in [nM]s
        // 1 / (1 + 1/(Math.pow(antibioticsCon / 100, 2)));
        double drugNowInNanoMolars = NgPerMlToNanomolars(antibioticsNow);
        double antibioticsEfficacy = 1 / ( 1 + (EC50 / StochasticDrug(drugNowInNanoMolars)));

        return antibioticsEfficacy;

    }

    double StochasticDrug(double drug){

        double stdDevOfGaussian = drug / 100;

        Rand random = new Rand();
        double stochasticDrug = random.Gaussian(drug, stdDevOfGaussian);
        stochasticDrug = stochasticDrug > 0.0 ? stochasticDrug : 0.0;

        return stochasticDrug;

    }

    double NgPerMlToNanomolars(double drugNow){
        return drugNow * Math.pow(10,3) / molarMassDrug;
    }

}

class NewExperiment extends AgentGrid2D<CartilageCells>{

    public int x = 100;
    public int y = 100;
    public int visScale = 3;
    public int numberOfTicksDelay;
    public int numberOfTicksDrug;
    public int numberOfTicks;
    public PDEGrid2D bacterialCon;
    public NeutrophilLayer neutrophilLayer;
    //public PDEGrid2D immuneResponseLevel; // similar to T-cell concentrations, but more generic
    public double antibioticsCon = 0;
    public double antibioticsConStomach = 0;
    public double corticosteroidCon = 0;
    public double corticosteroidConStomach = 0;
    public Rand rn;
    public double[] cellularBacterialCon = new double[length];
    public double[] cellularImmuneResponseLevel = new double[length];

    public double fixedDamageRate;

    public double stapyloReproductionRate = Math.pow(10,-4);
    public double damageRate = Math.pow(10,-6);
    public double deathProb = 7.02 * Math.pow(10,-4); // P_D
    public double staphyloDiffCoeff; // D_V [sigma^2 / min]

    AntibioticsDrug antibiotics;
    CorticosteroidDrug corticosteroid;

    String tabletsOrIntravenous = "tablets";

    public double immuneResponseDecay = 0.00005;
    public double immuneResponseDiffCoeff = 0.1;

    public boolean isAntibiotics = true;
    public boolean isCorticosteroid = true;

    public FileIO outFile;
    public FileIO paramFile;
    public FileIO concentrationsFile;
    public String outputDir;

    public NewExperiment(int xDim, int yDim, int visScale, Rand rn, boolean isAntibiotics, boolean isCorticosteroid, int numberOfTicksDelay, double staphyloDiffCoeff, String tabletsOrIntravenous, double fixedDamageRate){

        super(xDim, yDim, CartilageCells.class);
        this.x = xDim;
        this.y = yDim;
        this.visScale = visScale;
        this.numberOfTicksDelay = numberOfTicksDelay;
        this.staphyloDiffCoeff = staphyloDiffCoeff;
        this.rn = rn;

        this.tabletsOrIntravenous = tabletsOrIntravenous;
        this.isAntibiotics = isAntibiotics;
        this.isCorticosteroid = isCorticosteroid;

        this.fixedDamageRate = fixedDamageRate;

        if (tabletsOrIntravenous.equals("tablets")) {

            this.antibiotics = new AntibioticsDrug();
            this.corticosteroid = new CorticosteroidDrug();
            this.numberOfTicksDrug = 2 * 24 * 60; // we administer antibiotics for 5 days, i.e. 5*24*60 minutes

        } else {

            this.antibiotics = new AntibioticsDrug(5.0, isAntibiotics);
            this.corticosteroid = new CorticosteroidDrug(5.0, isCorticosteroid);
            this.numberOfTicksDrug = 4 * 24 * 60; // 4 days

        }

        bacterialCon = new PDEGrid2D(xDim, yDim);
        bacterialCon.Update();
        neutrophilLayer = new NeutrophilLayer(xDim, yDim);

        outputDir = this.OutputDirectory();
        outFile = new FileIO(outputDir.concat("/").concat("Out").concat(".csv"), "w");
        paramFile = new FileIO(outputDir.concat("/").concat("Param").concat(".csv"), "w");
        concentrationsFile = new FileIO(outputDir.concat("/").concat("concentrations").concat(".csv"), "w");

    }

    public void Init(){

        WriteHeader();
        String bacteriumSource = "injury"; // blood or injury

        if (bacteriumSource.equals("blood")) {
            int initialPlace1 = x / 2 - x % 2;
            int initialPlace2 = (99 * length) / 100 - (99 * length) % 100 + x / 2 - x % 2;
            for (int i = 0; i < length; i++) {
                CartilageCells c = NewAgentSQ(i);
                c.CellInit(true, false, false);
                if (i == initialPlace1 || i == initialPlace2) {
                    bacterialCon.Add(c.Isq(), 50.0);
                    bacterialCon.Update();
                }
            }
        }else if(bacteriumSource.equals("injury")){
            Rand random = new Rand();
            int initialPlace =  x / 2 - x % 2;
            for (int i = 0; i < length; i++) {
                CartilageCells c = NewAgentSQ(i);
                c.CellInit(true, false, false);
                if (i == initialPlace) {
                    bacterialCon.Add(c.Isq(), 50.0);
                    bacterialCon.Update();
                }
            }
        }
    }

    public double RunExperiment(GridWindow win, OpenGL2DWindow neutrophilWindow){

        double[] cellCounts = CountCells();
        // System.out.println(cellCounts[0]+", " + cellCounts[1] + ", " + cellCounts[2]);

        for (int tick = 0; tick < this.numberOfTicks; tick ++){

            if (neutrophilWindow.IsClosed()){
                break;
            }
            neutrophilWindow.TickPause(20);

            if ((numberOfTicksDelay == CorticosteroidsAntibiotics.BIG_VALUE) && (fixedDamageRate <= 100) && (((cellCounts[1] + cellCounts[2]) / 40000) * 100 >= fixedDamageRate)){
                this.numberOfTicksDelay = tick - 1;
                this.numberOfTicks = this.numberOfTicksDelay + this.numberOfTicksDrug;
                System.out.println("Diffusion coeff.: " + this.staphyloDiffCoeff + ". Damage info: " + fixedDamageRate + " percent damage was found at tick " + numberOfTicksDelay + ".");
            }

            TimeStep(tick);
            DrawModel(win);
            neutrophilLayer.DrawModel(neutrophilWindow);

            if( tick > 0 && ( (tick % (3*60)) == 0 ))
                win.ToPNG(outputDir + "3_hour_period_" + Integer.toString(tick/(3*60)) + ".jpg");

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
        for (CartilageCells cell: this){
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

    void TimeStepNeutrophils(int tick){

        // decay of the immuneResponseLevel
        double decayFactor = immuneResponseDecay + 0.01 * corticosteroid.CorticosteroidEfficacy(corticosteroidCon);
        for(Neutrophil cell: neutrophilLayer){
            cell.NeutrophilDecay(decayFactor);
        }

        // neutrophils arrive when bacterial concentration is high or when other neutrophils call them
        double corticosteroidSuppression = (1 - corticosteroid.CorticosteroidEfficacy(corticosteroidCon));
        double bacterialFactor = 1 / (1 + 1/(Math.pow(TotalBacterialCon(),2)));
        double arrivalFactor = bacterialFactor * corticosteroidSuppression;
        NeutrophilArrivalToSiteCalledByBacteria(arrivalFactor);

        // neutrophils "call" new ones by signaling, so this could be reasonable
        for(Neutrophil cell: neutrophilLayer){
            Rand rng = new Rand();
            if(rng.Double() < neutrophilLayer.signalingIntensity) {
                cell.NeutrophilArrivalToSiteCalledByNeutrophils(neutrophilLayer.signalSuccess * arrivalFactor);
            }
        }

        for(Neutrophil cell: neutrophilLayer){
            cell.NeutrophilMove();
        }
    }

    public void NeutrophilArrivalToSiteCalledByBacteria(double arrivalFactor){
        for (int i = 0; i < arrivalFactor; i++) {
            int[] arrivalPlace = Neutrophil.arrivalFromBorders();
            neutrophilLayer.NewAgentPT(arrivalPlace[0], arrivalPlace[1]).Init();
        }
    }
    void TimeStepStaphylo(int tick){

        // bacterial reproduction
        for (CartilageCells cell : this){
            bacterialCon.Add(cell.Isq(), stapyloReproductionRate * bacterialCon.Get(cell.Isq()));
        }

        bacterialCon.DiffusionADI(staphyloDiffCoeff);
        bacterialCon.Update();

        // removal of bacteria
        for (CartilageCells cell : this){
            // double removalEfficacy = 2/(1+Math.exp(100*drugNow));
            // double removalEfficacy = 100*Math.pow(drugNow, 2)/(1+100*Math.pow(drugNow,2));
            double drugBacterialRemovalEff = 0.01 * antibiotics.AntibioticsEfficacy(antibioticsCon);
            System.out.println(drugBacterialRemovalEff);
            double immuneBacterialRemovalEff = 0.5 * 1 / (1 + 1/(Math.pow(neutrophilLayer.PopAt(cell.Isq()), 2)));
            bacterialCon.Add(cell.Isq(), -drugBacterialRemovalEff * bacterialCon.Get(cell.Isq()));
            bacterialCon.Add(cell.Isq(), -immuneBacterialRemovalEff * bacterialCon.Get(cell.Isq()));
        }
        bacterialCon.Update();

    }

    void TimeStepDrug(int tick){

        if (this.tabletsOrIntravenous.equals("intravenous")){

            this.antibioticsCon = this.antibiotics.intravenousAntibioticsCon;
            this.corticosteroidCon = this.corticosteroid.intravenousCorticosteroidCon;

        } else if (this.tabletsOrIntravenous.equals("tablets")){

            // decay of the drug
            this.antibioticsCon -= this.antibiotics.drugDecay * this.antibioticsCon;
            this.corticosteroidCon -= this.corticosteroid.drugDecay * this.corticosteroidCon;

            // decay of the drug in the stomach
            // and appearance of the drug at the lung epithelial cells
            double transferQuantity = this.antibiotics.drugDecayStomach * this.antibioticsConStomach;
            this.antibioticsCon += transferQuantity;
            this.antibioticsConStomach -= transferQuantity;

            transferQuantity = this.corticosteroid.drugDecayStomach * this.corticosteroidConStomach;
            this.corticosteroidCon += transferQuantity;
            this.corticosteroidConStomach -= transferQuantity;

            // drug appearance in the stomach
            this.antibioticsConStomach += AntibioticsSourceStomach(tick);
            this.corticosteroidConStomach += CorticosteroidSourceStomach(tick);

        } else {
            System.out.println("Intravenous and tablets are the only two choices currently.");
        }

    }

    void TimeStepCartilageCells(int tick){

        for (CartilageCells cell : this){
            cell.CellDamage();
        }
        for (CartilageCells cell : this) {
            cell.CellDeath();
        }

    }

    void TimeStep(int tick){

        TimeStepStaphylo(tick);
        TimeStepCartilageCells(tick);
        TimeStepNeutrophils(tick); // immune response
        TimeStepDrug(tick); // both antibiotics and steroid

    }

    double TotalBacterialCon() {

        double totalBacterialCon = 0;
        for (int i = 0; i < length; i++){
            cellularBacterialCon[i] = bacterialCon.Get(i);
        }
        for (double bacterialConInCell : cellularBacterialCon ){
            totalBacterialCon = totalBacterialCon + bacterialConInCell;
        }
        return totalBacterialCon;
    }

    double AntibioticsSourceStomach(int tick){

        if ((tick > numberOfTicksDelay) && (isAntibiotics == true) && (((tick - numberOfTicksDelay) % (12 * 60)) == 1)) {
            return this.antibiotics.drugSourceStomach;
        } else {
            return 0.0;
        }
    }

    double CorticosteroidSourceStomach(int tick){

        if ((tick > numberOfTicksDelay) && (isCorticosteroid == true) && (((tick - numberOfTicksDelay) % (12 * 60)) == 1)) {
            return this.corticosteroid.drugSourceStomach;
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
        } else if (this.tabletsOrIntravenous.equals("intravenous")) {
            drugInfo = this.antibiotics.NgPerMlToNanomolars(this.antibiotics.intravenousAntibioticsCon);
        } else {
            drugInfo = this.antibiotics.drugSourceStomach;
        }

        String outputDir = projPath + "/" + date_time + this.tabletsOrIntravenous + drugInfo + "__diff" + this.staphyloDiffCoeff;
        if(this.fixedDamageRate < 100.0){
            outputDir +=  "__damagerate" + this.fixedDamageRate + "/";
        } else {
            outputDir +=  "__delaytime" + this.numberOfTicksDelay + "/";
        }

        new File(outputDir).mkdirs();

        return outputDir;

    }

    void DrawModel(GridWindow vis){

        for (int i = 0; i < length; i++) {
            CartilageCells drawMe = GetAgent(i);

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

            vis.SetPix(ItoX(i) + xDim, ItoY(i), HeatMapGRB(100*bacterialCon.Get(i)));

            vis.SetPix(ItoX(i) + 2 * xDim, ItoY(i), HeatMapBGR(0.1*neutrophilLayer.PopAt(i)));
        }
    }

}

class CartilageCells extends AgentSQ2Dunstackable<NewExperiment>{
    int CellType;

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

    public void CellDamage(){

        // nor antibiotics nor corticosteroids help the cells directly
        // difference w.r.t. certain antivirals, which prevent the cells from getting infected

        // chondral damage: depends on the immune response strength (a strong immune response kills bacteria, but damages the cartilage as well)

        double immuneConAtCell = G.neutrophilLayer.PopAt(Isq());
        double effectiveDamageProb = immuneConAtCell * G.damageRate * G.xDim * G.yDim;

        if (this.CellType == 0){ // healthy cell
            if (G.rn.Double() < effectiveDamageProb) {
                this.CellType = 1;
            }
        }
    }

    public void CellDeath(){

        if (this.CellType == 1) { // damaged
            if(G.rn.Double() < G.deathProb){
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
            return;
        }
    }

    public void NeutrophilArrivalToSiteCalledByNeutrophils(double signalSuccess){
        if((G.PopAt(Isq())<5) && G.rng.Double()<signalSuccess){
            int[] arrivalPlace = arrivalFromBorders();
            G.NewAgentPT(arrivalPlace[0], arrivalPlace[1]).Init();
        }
    }

    public void NeutrophilMove(){
        G.rng.RandomPointInCircle(0.5,G.moveCoords);
        MoveSafePT(Xpt()+G.moveCoords[0], Ypt()+G.moveCoords[1]);
    }

    public static int[] arrivalFromBorders(){
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



    //public void runNeutrophils(){

    //OpenGL2DWindow win = new OpenGL2DWindow("Neutrophils", 500, 500, x,y);

    // init model

    //for (int i = 0; i<timesteps; i++){
    //if (win.IsClosed()){
    //    break;
    //}
    //win.TickPause(20);
    //if(model.Pop()==0){
    //    model.NewAgentSQ(model.xDim/2, model.yDim/2).Init();
    //}

    //draw
    //model.DrawModel(win);
    //}
    //win.Close();

    //}
}