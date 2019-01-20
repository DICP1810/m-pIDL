/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package sixduplex;
import java.util.*;
import java.io.*;
/**
 *
 * @author shan
 */
//revised on 2012/12/03
//%(M + 1) = (1.11 x #C atoms) + (0.016 x #H atoms) + (0.38 x #N atoms) +0.04*o#+0.78*s#
//%(M + 2) @ [%(M + 1)2 / 200] + (0.20 x #O atoms) + (4.4 x #S atoms)
//http://www.chemistry.ccsu.edu/glagovich/teaching/316/ms/isotope.html

//what is the one unit of ms difference, for isotope? neutron,1.0087
//side-chain loss has not been considered
//how to calculate isotope content or ratios?
//the effect of isotope overlap for paired light and heavy fragment is not calculated yet
//the overlap between different light and heavy ion is considered and removed. can be refined durther, as long as the difference is not unit, still can be used?
//double ion threshold(M++>200?)
//thi is 为周愿设计,双甲醛标记
class NoDataColumnException extends Exception{}
public class RelativeQuant {
    List<HashMap<Character,Double>> allVarMod=new ArrayList<>();
    List<HashMap<Character,Double>> allFixMod=new ArrayList<>();
    List<PeptideFragmentation> allPep=new ArrayList<>();
    int multiple=0;
    double maxIonMass=200;

    double doubleIonThres=200;
    boolean removeIsotope=false;
    double mzThreshold=1000;
    //
    public ReadScanSummary expPeptide;
    //

    public List<List<String>> ionNames=new ArrayList<>(); 
    public List<List<List<Double>>> ionCounts=new ArrayList<>(); 
    public List<List<Double>> totalIonCounts=new ArrayList<>();
    public List<String> scans=new ArrayList<>();
     //
   
    File summaryFile;
    File mgfFile;
    //the precision of instrument for MSMS ion
    double threshold=0.01;
    //mass difference between light and heavy ions,should be integer unit or the mass of neutron?
    double setMassDiff=2.01255;//two neutron weight
    public RelativeQuant(File scansummaryfile,File mgffile){
        expPeptide=new ReadScanSummary(scansummaryfile);
        expPeptide.findExpMSMS(mgffile);
    }
    public void setThreshold(double thres,double massDiff){
        threshold=thres;
        setMassDiff=massDiff;
    }
    public void writeData(File scanSummaryFile){
        try{
            BufferedReader csvFile=new BufferedReader(new FileReader(scanSummaryFile));
            String dataRow;
            List<String> allData=new ArrayList<>();
            while ((dataRow=csvFile.readLine())!=null){
                allData.add(dataRow);
            }
            csvFile.close();
            try{
                FileWriter outFile = new FileWriter(scanSummaryFile);
                PrintWriter out = new PrintWriter(outFile);
                int line=totalIonCounts.size();
                String tem="";
                for(int i=0;i<multiple;i++)
                    tem=tem+"\t"+"exp"+Integer.toString(i+1);
                
                out.println(allData.get(0)+"\t"+"Ion Number"+"\t"+"Ion Name"+tem);
                
                for(int i=0;i<line;i++){
                    String te="";
                    for(int j=0;j<multiple;j++)
                        te=te+"\t"+Double.toString(totalIonCounts.get(i).get(j));
                    out.println(allData.get(i+1)+"\t"+Integer.toString(ionNames.get(i).size())+"\t"+ionNames.get(i)+te);
                }
                out.close();
            }
            catch (IOException e){
                e.printStackTrace();
            }
        }
        catch(FileNotFoundException e){
            System.out.println("file not found");
        }
        catch(IOException e){
       
            System.out.println("io exception");
        }
        
    }
        
            
    public void setAllMod(List<HashMap<Character,Double>> allvarmod, List<HashMap<Character,Double>> allfixmod){
        allVarMod=allvarmod;
        allFixMod=allfixmod;
        multiple=allVarMod.size();
        //System.out.println(multiple);
        for(int i=0;i<multiple;i++){
            PeptideFragmentation onePep=new PeptideFragmentation(allVarMod.get(i),allFixMod.get(i));
            //System.out.println(allVarMod.get(i));
            allPep.add(onePep);
        }
    }
    //transfer modified sequence from _DM(ox)K_ to 2.031.0
    public List<String> transferString(List<String> pepVarMod){
        int n=pepVarMod.size();
        //important: the modification must be in the same format
        List<String> pepNewVarMod=new ArrayList<>();
        /*for(int i=0;i<n;i++){
            String str1,str2,str3;
            str1=pepVarMod.get(i).replace("M(ox)", "3");
            if(str1.charAt(0)=='_'){
                str2=str1.replaceFirst("_", "2.");
            }    
            else
                str2="0."+str1;
            str3=str2.replace("K_","1")+".0";
            str3=str3.replace("K","1");
            pepNewVarMod.add(str3.replaceAll("[A-Z]", "0"));
        }*/
        //for six duplex modification
        
        for(int i=0;i<n;i++){
            String str1,str2,str3,str4;
            str1=pepVarMod.get(i).replace("M(ox)", "3");
            if(str1.charAt(0)=='_'){
                str2=str1.replace("_(di)", "1.");
                str2=str2.replace("_(ac)", "4.");
                str2=str2.replace("_","0.");
            }
            else
                str2="0."+str1;
            str3=str2.replace("K(di)_","2")+".0";
            str4=str3.replace("K(di)","2");
            //if end is not K
            str4=str4.replace("_","");
            //
            pepNewVarMod.add(str4.replaceAll("[A-Z]", "0"));
            
        }
        return pepNewVarMod;
    }
    public void calculateAllLightHeavy() throws NoDataColumnException{
        if(expPeptide.scanSummary.get("Modified sequence")==null||expPeptide.scanSummary.get("Sequence")==null||expPeptide.scanSummary.get("Scan number")==null||expPeptide.scanSummary.get("Charge")==null)
            throw new NoDataColumnException();
        List<String> pepVarMods=transferString(expPeptide.scanSummary.get("Modified sequence"));
        int n=pepVarMods.size();
        List<String> pepSeq=expPeptide.scanSummary.get("Sequence");
        List<String> pepScan=expPeptide.scanSummary.get("Scan number");
        List<String> pepCharge=expPeptide.scanSummary.get("Charge");
        
        String pepDFile;
        int charge;
        
        //some problems exists in reading charge because the previous column of it has ,, now change the column of charge to be after modified sequence in the csv file
        
        for(int i=0;i<n;i++){
            
            if(pepVarMods.get(i).length()>4){
                charge=Integer.parseInt(pepCharge.get(i));
                //to avoid big negative values
                if(charge<0)
                    charge=1;
                pepDFile=pepScan.get(i);
                //System.out.println(pepDFile);
                //get quantable light and heavy ion intensitied
                calculateLightHeavy(pepSeq.get(i),pepVarMods.get(i),charge,pepDFile,expPeptide.expMSMS.get(pepDFile));
                //System.out.println(i);
            }
        }
        
        System.out.println("success");
    }
    public boolean isAminoValid(String pep){
        String aminoAcidAbbreviation = "ARNDCEQGHILKMFPSTWYV";
        List<Character> amino=new ArrayList<>();
        for(int i=0;i<aminoAcidAbbreviation.length();i++)
            amino.add(aminoAcidAbbreviation.charAt(i));
        boolean isValid=true;
        for(int i=0;i<pep.length();i++){
            if(!amino.contains(pep.charAt(i))){
                isValid=false;
                break;
            }
        }
        return isValid;
            
    }
    //calculate theoretical ms/ms for light and heavy peptide,compare with expmsms,then determine which ion is quantifiable
    public void calculateLightHeavy(String peptidesequence, String peptidevarmod,int chargeprecursor,String dtafilename,HashMap<Double,Double> expmsms){
        List<String> ionName=new ArrayList<>();
        List<Double> allIonTotal=new ArrayList<>();
        boolean isValid=isAminoValid(peptidesequence);
        System.out.println(isValid);
        //judge the existence of valid amino acid
        if(isValid){
        for(int i=0;i<multiple;i++){
            allPep.get(i).calculateMSMS(peptidesequence, peptidevarmod,chargeprecursor);
            //System.out.println(allPep.get(i).varMod);
        }
        List<List<Double>> allTheoMz=new ArrayList<>();
        //theoretical ms/ms ion mz
        for(int i=0;i<multiple;i++){
            List<Double> oneTheoMz=new ArrayList<>(allPep.get(i).getTheoMSMS().values());
            allTheoMz.add(oneTheoMz);
            //System.out.println(oneTheoMz);
        }
        //theoretical ms/ms ion name
        List<String> theoName=new ArrayList<>(allPep.get(0).getTheoMSMS().keySet());
        System.out.println(theoName);
        List<List<Double>> ionCount=new ArrayList<>();
        for(int i=0;i<multiple;i++){
            List<Double> oneIonCount=new ArrayList<>();
            ionCount.add(oneIonCount);
        }
        List<Integer> findIndex=new ArrayList<>(findCommonIon(ionCount,allTheoMz, expmsms,dtafilename));
        int num=findIndex.size();

        if(num>0){
            List<List<Double>> newAllTheoMz=new ArrayList<>();
            for(int i=0;i<multiple;i++){
                List<Double> newOneTheoMz=new ArrayList<>();
                newAllTheoMz.add(newOneTheoMz);
            }
            List<String> newNames=new ArrayList<>();
            List<Integer> newCharge=new ArrayList<>();
            for(int i=0;i<findIndex.size();i++){
                for(int j=0;j<multiple;j++)
                    newAllTheoMz.get(j).add(allTheoMz.get(j).get(findIndex.get(i)));
                newNames.add(theoName.get(findIndex.get(i)));
                newCharge.add(allPep.get(0).getChargeMS().get(findIndex.get(i)));
            }
            
            List<Double> ionMass=new ArrayList<>();
            List<List<Double>> quantAllIonCount=new ArrayList<>();
            for(int i=0;i<multiple;i++){
                List<Double> quantOneIonCount=new ArrayList<>();
                quantAllIonCount.add(quantOneIonCount);
            }
            
            for(int i=0;i<multiple;i++)
                allIonTotal.add(0.0);
            for(int i=0;i<findIndex.size();i++){
                ionName.add(newNames.get(i));
                ionMass.add(newAllTheoMz.get(0).get(i));
                for(int j=0;j<multiple;j++){
                    //each ions for one spectrum of one labeling
                    quantAllIonCount.get(j).add(ionCount.get(j).get(i));
                    //total ions for one spectrum of one labeling
                    allIonTotal.set(j,allIonTotal.get(j)+ionCount.get(j).get(i));
                }
            }
            ionNames.add(ionName);
            scans.add(dtafilename);
            totalIonCounts.add(allIonTotal);
            ionCounts.add(quantAllIonCount);
        }
        else{
            ionNames.add(ionName);
            for(int i=0;i<multiple;i++)
                allIonTotal.add(0.0);
            totalIonCounts.add(allIonTotal);
        }
        }
        else{
            ionNames.add(ionName);
            for(int i=0;i<multiple;i++)
                allIonTotal.add(0.0);
            totalIonCounts.add(allIonTotal);
        }
        
    }
     //find the common ion intensity for different labeling
    public List<Integer> findCommonIon(List<List<Double>> ionCount,List<List<Double>> alltheomz, HashMap<Double,Double> expmsms,String dtafilename){
        int charge;
        int expL=expmsms.size();
        int theoL=alltheomz.get(0).size();
        double temp;
        List<Integer> findIndex=new ArrayList<>();
        List<Double> expMass=new ArrayList<>(expmsms.keySet());
        
        for(int i=0;i<multiple;i++)
            System.out.println(alltheomz.get(i).size());
        for(int i=0;i<theoL;i++){
                List<Double> commonIonIntensity=new ArrayList<>();
                for(int j=0;j<multiple;j++)
                    commonIonIntensity.add(0.0);
                charge=allPep.get(0).getChargeMS().get(i);
                boolean flag=false;
                for(int j=0;j<multiple;j++){
                    if(alltheomz.get(j).get(i)<doubleIonThres){
                        flag=true;
                        break;
                    }
                }
                if(charge==2&&flag)
                    ;
                else{
                    for(int j=0;j<expL;j++){
                        temp=expMass.get(j);
                        //temporary to limit only the small ions
                        
                        for(int k=0;k<multiple;k++){
                            //System.out.println(alltheomz.get(k).get(i));
                            if(Math.abs(temp-alltheomz.get(k).get(i))<threshold&&alltheomz.get(k).get(i)<mzThreshold){
                                //System.out.println(alltheomz.get(k).get(i));
                                //System.out.println(temp);
                                //System.out.println(expmsms.get(temp));
                                if(expmsms.get(temp)>commonIonIntensity.get(k)&& alltheomz.get(k).get(i)<maxIonMass)
                                    commonIonIntensity.set(k,expmsms.get(temp));
                            }
                        }
                    }
                    //System.out.println(i);
                    boolean flag1=true;
                    for(int j=0;j<multiple;j++){
                        if(commonIonIntensity.get(j)<=0){
                            flag1=false;
                            break;
                        }
                    }
                     //System.out.println(commonIonIntensity);
                     //System.out.println(flag1);
                    if(flag1){
                        findIndex.add(i);
                        System.out.println(commonIonIntensity);
                        for(int j=0;j<multiple;j++){
                            ionCount.get(j).add(commonIonIntensity.get(j));
                        }
                    }
                }
        }
        return findIndex;
    }
      
    
}
