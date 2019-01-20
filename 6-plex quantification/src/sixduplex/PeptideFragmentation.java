/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package sixduplex;

import java.util.*;
//compare with http://db.systemsbiology.net:8080/proteomicsToolkit/FragIonServlet.html
//http://www.matrixscience.com/help/fragmentation_help.html
//according to literature, only list neutral loss of h2o and nh3 for a,b,y ions.
//x ion is not listed. nth fragment ion is listed for a,b and z(details see in the program).
// the loss of side chain has not been considered
//double ion occurs if molecular weight exceeds certain limit, here is set as 400
//http://www.ionsource.com/Card/aatable/aatable.htm: amino acid weight
//proton molweight:1.0073.
/*
 * @author shan
 */
//revised on 2012/12/03
class NoAminoAcidException extends Exception{}
public class PeptideFragmentation {
    
    ////to control the mass of the theoretical mass to be used
    double maxIonMass=10000;
    ////to judge the last amino of pep
    boolean isK=false;
    ///
    double hydrogenMolWeight=1.0078250;
    double protonMolWeight=1.0073;
    double oxygenMolWeight=15.994915;
    double waterMW=2*hydrogenMolWeight+oxygenMolWeight;
    double carbonMolWeight=12.0000;
    double nitrogenMolWeight=14.003074;
    double ammoniaMW=3*hydrogenMolWeight+nitrogenMolWeight;
    double coMW=carbonMolWeight+oxygenMolWeight;
    double[]aminoAcidResMolWeight={71.037114,156.101111,114.042927,115.026943,103.009185,129.042593,128.058578,57.021464,137.058912,113.084064,113.084064,128.094963,131.040485,147.068414,97.052764,87.032028,101.047679,186.079313,163.06332,99.068414};
    HashMap<Character,Double> aminoAcidResidueMolWeight=new HashMap<>();
    String aminoAcidAbbreviation = "ARNDCEQGHILKMFPSTWYV";
    
    private boolean aIonSwitch=true;
    private boolean xIonSwitch=false;
    private boolean bIonSwitch=true;
    private boolean yIonSwitch=true;
    private boolean cIonSwitch=false;
    private boolean zIonSwitch=false;
    private boolean lossWaterSwitch=false;
    private boolean lossAmmoniaSwitch=false;
    
    private boolean doubleFragIon=false;
    private boolean motherIonSwitch=false;
    private boolean nthIonSwitch=false;        
    public HashMap<String,Double> theoMSMS=null;
    
    public List<Integer> chargeMS=new ArrayList<>();
    //double[] bIon;
    List<Double> bIon=new ArrayList<>();
    
    int lengthPeptide;
    int chargePrecursor;
    public double peptideMass;
    String peptideSequence;
    String peptideVarMod;
    
    HashMap<Character,Double> varMod=new LinkedHashMap<>();
    HashMap<Character,Double> fixedMod=new LinkedHashMap<>();;
    
    public PeptideFragmentation(HashMap<Character,Double> varmod,HashMap<Character,Double> fixedmod){
        varMod=varmod;
        fixedMod=fixedmod;
        theoMSMS=new LinkedHashMap<>();
    }
    public void calculateMSMS(String peptidesequence,String peptidevarmod,int chargeprecursor){
        peptideInformation(peptidesequence,peptidevarmod,chargeprecursor);
        try{
            calculateBIon();
            calculatePepMass();
            //theoMSMS.clear();
            //chargeMS.clear();
            theoMSMS=new LinkedHashMap<>();
            chargeMS=new ArrayList<>();
            calculateMotherIon();
            calculateOtherIon();
        }
        catch(NoAminoAcidException e){
            System.out.println("amino acid in "+peptidesequence+" doesn't exist");
            e.printStackTrace(System.out);
        }
    }
    public void peptideInformation(String peptidesequence,String peptidevarmod,int chargeprecursor){
        peptideSequence=peptidesequence;
        if(peptideSequence.charAt(peptideSequence.length()-1)=='K')
            isK=true;
        else
            isK=false;
        peptideVarMod=peptidevarmod;
        chargePrecursor=chargeprecursor;
        lengthPeptide=peptideSequence.length();
    }
    public void setPeptideMass(double peptidemass){
        peptideMass=peptidemass;
    }
    public void setSwitch(HashMap<String,Boolean> ionSwitch){
        aIonSwitch=ionSwitch.get("a");
        xIonSwitch=ionSwitch.get("x");
        bIonSwitch=ionSwitch.get("b");
        yIonSwitch=ionSwitch.get("y");
        cIonSwitch=ionSwitch.get("c");
        zIonSwitch=ionSwitch.get("z");
        lossWaterSwitch=ionSwitch.get("h2o");
        lossAmmoniaSwitch=ionSwitch.get("nh3");
        doubleFragIon=ionSwitch.get("doubleIon");
        motherIonSwitch=ionSwitch.get("motherIon");
        nthIonSwitch=ionSwitch.get("nthIon");
    }        
    public void calculateMotherIon(){
        if(motherIonSwitch){
        theoMSMS.put("M+",peptideMass+protonMolWeight);
        chargeMS.add(1);
        if(chargePrecursor>=2){ 
            theoMSMS.put("M++",(peptideMass+2.0*protonMolWeight)/2.0);
            chargeMS.add(2);
        }
        if(chargePrecursor==3){
            theoMSMS.put("M++",(peptideMass+3.0*protonMolWeight)/3.0);
            chargeMS.add(3);
        }
        
        if(lossWaterSwitch){
            theoMSMS.put("M+LossWater",peptideMass+protonMolWeight-waterMW);
            chargeMS.add(1);
            if(chargePrecursor>=2){ 
                theoMSMS.put("M++LossWater",(peptideMass+2.0*protonMolWeight-waterMW)/2.0);
                chargeMS.add(2);
            }
            if(chargePrecursor==3){
                theoMSMS.put("M+++LossWater",(peptideMass+3.0*protonMolWeight-waterMW)/3.0);
                chargeMS.add(3);
            }
        }
        if(lossAmmoniaSwitch){
            theoMSMS.put("M+LossAmmonia",peptideMass+protonMolWeight-ammoniaMW);
            chargeMS.add(1);
            if(chargePrecursor>=2){ 
                theoMSMS.put("M++LossAmmonia",(peptideMass+2.0*protonMolWeight-ammoniaMW)/2.0);
                chargeMS.add(2);
            }
            if(chargePrecursor==3){
                theoMSMS.put("M+++LossAmmonia",(peptideMass+3.0*protonMolWeight-ammoniaMW)/3.0);
                chargeMS.add(3);
            }
        }
        }    
    }
            
    public void calculateBIon()throws NoAminoAcidException{
        bIon=new ArrayList<>();
        //below calculate b ion molecular weight
        //substring, the first param is the start index, the last param is the last index+1
        System.out.println(peptideSequence);
        //System.out.println(peptideVarMod);
        
        int len=peptideVarMod.length();
        double nTerminalModMass=0;
        if(len>4&&varMod.containsKey(peptideVarMod.charAt(0)))
            nTerminalModMass+=varMod.get(peptideVarMod.charAt(0));
        
        //calculate MolWeight of amino acid residue after fixed modification
        int index=0;
        for(char eachAminoAcid: aminoAcidAbbreviation.toCharArray())
        {
            //residue mass equals to amino acid mass minus water mol weight
            if(fixedMod.containsKey(eachAminoAcid))
                //if fixed modification exists, the molecular weight of the amino acid residue should be adjusted
                aminoAcidResidueMolWeight.put(eachAminoAcid,aminoAcidResMolWeight[index]+fixedMod.get(eachAminoAcid));
    
            else
                aminoAcidResidueMolWeight.put(eachAminoAcid,aminoAcidResMolWeight[index]);
            
            index+=1;
        }
        String peptideVariableModificationNonTerminal="";
        if(len>4){
            peptideVariableModificationNonTerminal=peptideVarMod.substring(2, peptideVarMod.length()-2);
            //System.out.println( peptideVariableModificationNonTerminal);
        }
        index=0;
        double temp;
        for(char eachAminoAcid: (peptideSequence.substring(0,peptideSequence.length()-1)).toCharArray())
        {
            if(aminoAcidResidueMolWeight.get(eachAminoAcid)==null)
                throw new NoAminoAcidException();
            else{
            if(index==0)
            {                
                //for the first b ion,the modification to the n-terminal should be considered
                bIon.add(aminoAcidResidueMolWeight.get(eachAminoAcid)+nTerminalModMass);
            }
            else
            {
                //System.out.println(bIon.get(index-1));
                //System.out.println(aminoAcidResidueMolWeight.get(eachAminoAcid));
                temp=bIon.get(index-1)+aminoAcidResidueMolWeight.get(eachAminoAcid);
                //if(index==21)
                        //System.out.println(peptideSequence);
                if(len>4&&varMod.containsKey(peptideVariableModificationNonTerminal.charAt(index))){
                    
                    temp=temp+varMod.get(peptideVariableModificationNonTerminal.charAt(index));
                }
                bIon.add(temp);
            }
            }
            index+=1;
        }
    }
    public void calculatePepMass(){
        double cTerminalModMass=0;
        //make sure there is variable modification or not by using len>4 to judge
        int len=peptideVarMod.length();
        if(len>4&&varMod.containsKey(peptideVarMod.charAt(peptideVarMod.length()-1)))
            cTerminalModMass=varMod.get(peptideVarMod.charAt(peptideVarMod.length()-1));
        double lastAminoAcidModMass=0;
        if(len>4&&varMod.containsKey(peptideVarMod.charAt(peptideVarMod.length()-3)))
            lastAminoAcidModMass=varMod.get(peptideVarMod.charAt(peptideVarMod.length()-3));
        //calculate the peptide mw 
        peptideMass=bIon.get(bIon.size()-1)+aminoAcidResidueMolWeight.get(peptideSequence.charAt(peptideSequence.length()-1))+lastAminoAcidModMass+cTerminalModMass+waterMW;    
        System.out.println(peptideMass);
    }
            
    public void calculateOtherIon()
    {
        //this option should be put after both light and heavy peptide ion has been fragmented.
        double doubleIonThres=0;
        int length=lengthPeptide-1;
        double temp;
        //for bn ion
        if(nthIonSwitch&&bIonSwitch){
            temp=peptideMass-waterMW+protonMolWeight;
            if(!lossWaterSwitch){ //if allow loss water bn equals to peptide loss water
                theoMSMS.put("b"+Integer.toString(lengthPeptide)+"+",temp);
                chargeMS.add(1);
            }
            if(lossWaterSwitch){
                theoMSMS.put("b"+Integer.toString(lengthPeptide)+"+"+"LossWater",temp-waterMW);
                chargeMS.add(1);
            }
            if(lossAmmoniaSwitch){
                theoMSMS.put("b"+Integer.toString(lengthPeptide)+"+"+"LossAmmonia",temp-ammoniaMW);
                chargeMS.add(1);
            }
            if(chargePrecursor>=2&&doubleFragIon){
                temp=peptideMass-waterMW+2.0*protonMolWeight;
                if(!lossWaterSwitch&&temp>doubleIonThres){
                    theoMSMS.put("b"+Integer.toString(lengthPeptide)+"++",temp/2.0);
                    chargeMS.add(2);
                }
                if(lossWaterSwitch&&temp-waterMW>doubleIonThres){
                    theoMSMS.put("b"+Integer.toString(lengthPeptide)+"++"+"LossWater",(temp-waterMW)/2.0);
                    chargeMS.add(2);
                }
                if(lossAmmoniaSwitch&&temp-ammoniaMW>doubleIonThres){
                    theoMSMS.put("b"+Integer.toString(lengthPeptide)+"++"+"LossAmmonia",(temp-ammoniaMW)/2.0);
                    chargeMS.add(2);
                }
            }
        }
        //bn ion end
        if(bIonSwitch){
            for(int i=0;i<length;i++){
                temp=bIon.get(i)+protonMolWeight;
                //limit only small mass ions
                if(temp<maxIonMass){
                  theoMSMS.put("b"+Integer.toString(i+1)+"+",temp);
                  chargeMS.add(1);
                }
                if(lossWaterSwitch){
                    theoMSMS.put("b"+Integer.toString(i+1)+"+"+"LossWater",temp-waterMW);
                    chargeMS.add(1);
                }
                if(lossAmmoniaSwitch){
                    theoMSMS.put("b"+Integer.toString(i+1)+"+"+"LossAmmonia",temp-ammoniaMW);
                    chargeMS.add(1);
                }
            }
            if(chargePrecursor>=2&&doubleFragIon){
                for(int i=0;i<length;i++){
                    temp=bIon.get(i)+2.0*protonMolWeight;
                    if(temp>doubleIonThres){
                        theoMSMS.put("b"+Integer.toString(i+1)+"++",temp/2.0);
                        chargeMS.add(2);
                    }
                    if(lossWaterSwitch&&temp-waterMW>doubleIonThres){
                        theoMSMS.put("b"+Integer.toString(i+1)+"++"+"LossWater",(temp-waterMW)/2.0);
                        chargeMS.add(2);
                    }
                    if(lossAmmoniaSwitch&&temp-ammoniaMW>doubleIonThres){
                        theoMSMS.put("b"+Integer.toString(i+1)+"++"+"LossAmmonia",(temp-ammoniaMW)/2.0);
                        chargeMS.add(2);
                    }
                }
            }
        }
        //yn ion equals to peptide mother       
        //
        if(!motherIonSwitch&&yIonSwitch&&nthIonSwitch){
        theoMSMS.put("M+",peptideMass+protonMolWeight);
        chargeMS.add(1);
        if(chargePrecursor>=2){ 
            theoMSMS.put("M++",(peptideMass+2.0*protonMolWeight)/2.0);
            chargeMS.add(2);
        }
        
        if(lossWaterSwitch){
            theoMSMS.put("M+LossWater",peptideMass+protonMolWeight-waterMW);
            chargeMS.add(1);
            if(chargePrecursor>=2){ 
                theoMSMS.put("M++LossWater",(peptideMass+2.0*protonMolWeight-waterMW)/2.0);
                chargeMS.add(2);
            }
        }
        if(lossAmmoniaSwitch){
            theoMSMS.put("M+LossAmmonia",peptideMass+protonMolWeight-ammoniaMW);
            chargeMS.add(1);
            if(chargePrecursor>=2){ 
                theoMSMS.put("M++LossAmmonia",(peptideMass+2.0*protonMolWeight-ammoniaMW)/2.0);
                chargeMS.add(2);
            }
        }
        }    

        if(yIonSwitch){
            //exclude y1 ions
            for(int i=1;i<length;i++){
                temp=peptideMass-bIon.get(length-i-1)+protonMolWeight;
                //limit small ions and only K-end peptide
                if(isK==true&&temp<maxIonMass){
                  theoMSMS.put("y"+Integer.toString(i+1)+"+",temp);
                  chargeMS.add(1);
                }
                if(lossWaterSwitch){
                    theoMSMS.put("y"+Integer.toString(i+1)+"+"+"LossWater",temp-waterMW);
                    chargeMS.add(1);
                }
                if(lossAmmoniaSwitch){
                    theoMSMS.put("y"+Integer.toString(i+1)+"+"+"LossAmmonia",temp-ammoniaMW);
                    chargeMS.add(1);
                }
            }
            if(chargePrecursor>=2&&doubleFragIon){
                for(int i=0;i<length;i++){
                    temp=peptideMass-bIon.get(length-i-1)+2.0*protonMolWeight;
                    if(temp>doubleIonThres){
                        theoMSMS.put("y"+Integer.toString(i+1)+"++",temp/2.0);
                        chargeMS.add(2);
                    }
                    if(lossWaterSwitch&&temp-waterMW>doubleIonThres){
                        theoMSMS.put("y"+Integer.toString(i+1)+"++"+"LossWater",(temp-waterMW)/2.0);
                        chargeMS.add(2);
                    }
                    if(lossAmmoniaSwitch&&temp-ammoniaMW>doubleIonThres){
                        theoMSMS.put("y"+Integer.toString(i+1)+"++"+"LossAmmonia",(temp-ammoniaMW)/2.0);
                        chargeMS.add(2);
                    }
                }
            }
                
        }
        //for An ion
        if(aIonSwitch&&nthIonSwitch){
            temp=peptideMass-coMW-waterMW+protonMolWeight;
            theoMSMS.put("a"+Integer.toString(lengthPeptide)+"+",temp);
            chargeMS.add(1);
            if(lossWaterSwitch){
                theoMSMS.put("a"+Integer.toString(lengthPeptide)+"+"+"LossWater",temp-waterMW);
                chargeMS.add(1);
            }
            if(lossAmmoniaSwitch){
                theoMSMS.put("a"+Integer.toString(lengthPeptide)+"+"+"LossAmmonia",temp-ammoniaMW);
                chargeMS.add(1);
            }
            if(chargePrecursor>=2&&doubleFragIon){
                temp=peptideMass-coMW-waterMW+2.0*protonMolWeight;
                if(temp>doubleIonThres){
                    theoMSMS.put("a"+Integer.toString(lengthPeptide)+"++",temp/2.0);
                    chargeMS.add(2);
                }
                if(lossWaterSwitch&&temp-waterMW>doubleIonThres){
                    theoMSMS.put("a"+Integer.toString(lengthPeptide)+"++"+"LossWater",(temp-waterMW)/2.0);
                    chargeMS.add(2);
                }
                if(lossAmmoniaSwitch&&temp-ammoniaMW>doubleIonThres){
                    theoMSMS.put("a"+Integer.toString(lengthPeptide)+"++"+"LossAmmonia",(temp-ammoniaMW)/2.0);
                    chargeMS.add(2);
                }
            }
        }
        //An end
        if(aIonSwitch){
            for(int i=0;i<length;i++){
                temp=bIon.get(i)+protonMolWeight-coMW;
                //limit only small ion
                if(temp<maxIonMass)
                {
                  theoMSMS.put("a"+Integer.toString(i+1)+"+",temp);
                  chargeMS.add(1);
                }
                if(lossWaterSwitch){
                    theoMSMS.put("a"+Integer.toString(i+1)+"+"+"LossWater",temp-waterMW);
                    chargeMS.add(1);
                }
                if(lossAmmoniaSwitch){
                    theoMSMS.put("a"+Integer.toString(i+1)+"+"+"LossAmmonia",temp-ammoniaMW);
                    chargeMS.add(1);
                }
            }
            if(chargePrecursor>=2&&doubleFragIon){
                for(int i=0;i<length;i++){
                    temp=bIon.get(i)+2.0*protonMolWeight-coMW;
                    if(temp>doubleIonThres){
                        theoMSMS.put("a"+Integer.toString(i+1)+"++",temp/2.0);
                        chargeMS.add(2);
                    }
                    if(lossWaterSwitch&&temp-waterMW>doubleIonThres){
                        theoMSMS.put("a"+Integer.toString(i+1)+"++"+"LossWater",(temp-waterMW)/2.0);
                        chargeMS.add(2);
                    }
                    if(lossAmmoniaSwitch&&temp-ammoniaMW>doubleIonThres){
                        theoMSMS.put("a"+Integer.toString(i+1)+"++"+"LossAmmonia",(temp-ammoniaMW)/2.0);
                        chargeMS.add(2);
                    }
                }
            }
        }
        //x ion doen't exist usually, so the switch is turned off
        if(xIonSwitch){
            for(int i=0;i<length;i++){
                temp=peptideMass-bIon.get(length-i-1)+protonMolWeight+coMW;
                theoMSMS.put("x"+Integer.toString(i+1)+"+",temp);
                chargeMS.add(1);
                if(lossWaterSwitch){
                    theoMSMS.put("x"+Integer.toString(i+1)+"+"+"LossWater",temp-waterMW);
                    chargeMS.add(1);
                }
                if(lossAmmoniaSwitch){
                    theoMSMS.put("x"+Integer.toString(i+1)+"+"+"LossAmmonia",temp-ammoniaMW);
                    chargeMS.add(1);
                }
            }
            if(chargePrecursor>=2&&doubleFragIon){
                for(int i=0;i<length;i++){
                    temp=peptideMass-bIon.get(length-i-1)+2.0*protonMolWeight+coMW;
                    if(temp>doubleIonThres){
                        theoMSMS.put("x"+Integer.toString(i+1)+"++",temp/2.0);
                        chargeMS.add(2);
                    }
                    if(lossWaterSwitch&&temp-waterMW>doubleIonThres){
                        theoMSMS.put("x"+Integer.toString(i+1)+"++"+"LossWater",(temp-waterMW)/2.0);
                        chargeMS.add(2);
                    }
                    if(lossAmmoniaSwitch&&temp-ammoniaMW>doubleIonThres){
                        theoMSMS.put("x"+Integer.toString(i+1)+"++"+"LossAmmonia",(temp-ammoniaMW)/2.0);
                        chargeMS.add(2);
                    }
                }
            }
                
        }
        //cn ion dosenot exist and is omitted
        if(cIonSwitch){
            for(int i=0;i<length;i++){
                temp=bIon.get(i)+protonMolWeight+ammoniaMW;
                theoMSMS.put("c"+Integer.toString(i+1)+"+",temp);
                chargeMS.add(1);
                /*
                if(lossWaterSwitch){
                    theoMSMS.put("c"+Integer.toString(i+1)+"+"+"LossWater",temp-waterMW);
                    chargeMS.add(1);
                }
                if(lossAmmoniaSwitch){
                    theoMSMS.put("c"+Integer.toString(i+1)+"+"+"LossAmmonia",temp-ammoniaMW);
                    chargeMS.add(1);
                }
                */
            }
            if(chargePrecursor>=2&&doubleFragIon){
                for(int i=0;i<length;i++){
                    temp=bIon.get(i)+2.0*protonMolWeight+ammoniaMW;
                    if(temp>doubleIonThres){
                        theoMSMS.put("c"+Integer.toString(i+1)+"++",temp/2.0);
                        chargeMS.add(2);
                    }
                    /*
                    if(lossWaterSwitch){
                        theoMSMS.put("c"+Integer.toString(i+1)+"++"+"LossWater",(temp-waterMW)/2.0);
                        chargeMS.add(2);
                    }
                    if(lossAmmoniaSwitch){
                        theoMSMS.put("c"+Integer.toString(i+1)+"++"+"LossAmmonia",(temp-ammoniaMW)/2.0);
                        chargeMS.add(2);
                    }
                    */
                }
            }
        }
        //zn begin
        if(zIonSwitch&&nthIonSwitch){
            temp=peptideMass+protonMolWeight-ammoniaMW;
            //if allow loss ammonia, zn equals to mother ion loss ammonia
            if(!lossAmmoniaSwitch){
                theoMSMS.put("z"+Integer.toString(lengthPeptide)+"+",temp);
                chargeMS.add(1);
            }
            /*
            if(lossWaterSwitch){
                theoMSMS.put("z"+Integer.toString(lengthPeptide)+"+"+"LossWater",temp-waterMW);
                chargeMS.add(1);
            }
            if(lossAmmoniaSwitch){
                theoMSMS.put("z"+Integer.toString(lengthPeptide)+"+"+"LossAmmonia",temp-ammoniaMW);
                chargeMS.add(1);
            }
            */
            if( chargePrecursor>=2&&doubleFragIon){
                temp=peptideMass+2.0*protonMolWeight-ammoniaMW;
                if(!lossAmmoniaSwitch&&temp>doubleIonThres){
                    theoMSMS.put("z"+Integer.toString(lengthPeptide)+"++",temp/2.0);
                    chargeMS.add(2);
                }
              /*
                if(lossWaterSwitch){
                    theoMSMS.put("z"+Integer.toString(lengthPeptide)+"++"+"LossWater",(temp-waterMW)/2.0);
                    chargeMS.add(2);
                }
                if(lossAmmoniaSwitch){
                    theoMSMS.put("z"+Integer.toString(lengthPeptide)+"++"+"LossAmmonia",(temp-ammoniaMW)/2.0);
                    chargeMS.add(2);
                }
                */
            }
                
        }
        //zn ion end
        //http://www.matrixscience.com/help/fragmentation_help.html
        //according to above link, c,z ion normally has not neutral loss
        if(zIonSwitch){
            for(int i=0;i<length;i++){
                temp=peptideMass-bIon.get(length-i-1)+protonMolWeight-ammoniaMW;
                if(!lossAmmoniaSwitch||!yIonSwitch){
                    theoMSMS.put("z"+Integer.toString(i+1)+"+",temp);
                    chargeMS.add(1);
                }
                /*
                if(lossWaterSwitch){
                    theoMSMS.put("z"+Integer.toString(i+1)+"+"+"LossWater",temp-waterMW);
                    chargeMS.add(1);
                }
                if(lossAmmoniaSwitch){
                    theoMSMS.put("z"+Integer.toString(i+1)+"+"+"LossAmmonia",temp-ammoniaMW);
                    chargeMS.add(1);
                }
                */
            }
            if(chargePrecursor>=2&&doubleFragIon){
                for(int i=0;i<length;i++){
                    temp=peptideMass-bIon.get(length-i-1)+2.0*protonMolWeight-ammoniaMW;
                    if(!lossAmmoniaSwitch||!yIonSwitch){
                        if(temp>doubleIonThres){
                            theoMSMS.put("z"+Integer.toString(i+1)+"++",temp/2.0);
                            chargeMS.add(2);
                        }
                    }
                    /*
                    if(lossWaterSwitch){
                        theoMSMS.put("z"+Integer.toString(i+1)+"++"+"LossWater",(temp-waterMW)/2.0);
                        chargeMS.add(2);
                    }
                    if(lossAmmoniaSwitch){
                        theoMSMS.put("z"+Integer.toString(i+1)+"++"+"LossAmmonia",(temp-ammoniaMW)/2.0);
                        chargeMS.add(2);
                    }
                    */
                }
            }
                
        }
    }
    public HashMap<String,Double> getTheoMSMS(){
        return theoMSMS;
    }
    public List<Integer> getChargeMS(){
        return chargeMS;
    }
    //donnot sort the theoretical msms values, otherwise the light and heavy ion is not corresponding.
    /*public LinkedHashMap sortHashMapByValues(HashMap<String,Double> passedMap) {
        List<String> mapKeys = new ArrayList(passedMap.keySet());
        List<Double> mapValues = new ArrayList(passedMap.values());
        Collections.sort(mapValues);
        Collections.sort(mapKeys);
        LinkedHashMap sortedMap = new LinkedHashMap();
        Iterator valueIt = mapValues.iterator();
        while (valueIt.hasNext()) {
            Object val = valueIt.next();
            Iterator keyIt = mapKeys.iterator();
            while (keyIt.hasNext()) {
                Object key = keyIt.next();
                String comp1 = passedMap.get(key).toString();
                String comp2 = val.toString();
                if (comp1.equals(comp2)){
                    passedMap.remove(key);
                    mapKeys.remove(key);
                    sortedMap.put((String)key, (Double)val);
                    break;
                }
            }
        }
        return sortedMap;
    }*/
}
//sometimes there is no modification, should be considered in the mass calculation
 
