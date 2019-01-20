/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package sixduplex;
import java.io.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 *
 * @author shan
 */
//revised on 2012/12/03
//read the scan summary file, find the exerimental spectrum for each scan
//due to the existing of non-unique peptides, in scansummary file,there is redundant dta file name, which doesn't affect the result
////需要定制的部分:读取scan summary文件时,需要读入扫描数(或者每个扫描的文件名),肽的序列,肽的修饰.首先是获得每列的标题,然后读入每行数据,根据扫描数
//读入mgf文件中队应的二级质谱图.

public class ReadScanSummary {
    //the data of scan summary file, containing all of the columns
    public HashMap<String,List<String>> scanSummary=new LinkedHashMap<>();
    private File scanSummaryFile;
    private File mgfFile;
    private HashMap<String,Integer> expMSMSBegin=new LinkedHashMap<>();
    private HashMap<String,Integer> expMSMSEnd=new LinkedHashMap<>();
    //store all of the experimental ms/ms listed in the scan summary file
    public HashMap<String,HashMap<Double,Double>> expMSMS=new LinkedHashMap<>();
    public HashMap<String,Integer> scanOrder=new LinkedHashMap<>();
    public  ReadScanSummary(File scansummaryfile){
        scanSummaryFile=scansummaryfile;
        readSummaryFile();
    }
    public void readSummaryFile(){
        scanSummary.clear();
        try{
            BufferedReader csvFile=new BufferedReader(new FileReader(scanSummaryFile));
            String dataRow = csvFile.readLine();
            //revised because the first line of scan summary file contains the title of each column 
            String[] dataArrayIndex=dataRow.split("\\t");
            //find the number of columns in scan summary file
            int n=dataArrayIndex.length;
            List<String> allData=new ArrayList<>();
            while ((dataRow=csvFile.readLine())!=null){
                allData.add(dataRow);
            }
            csvFile.close();
            String[] dataArray;
            String temp;
            int len=allData.get(0).split("\\t").length;
            if(len<n)
                n=len;
            for (int i=0;i<n;i++){
                    List<String> eachColumn=new ArrayList<>();
                    for(int j=0;j<allData.size();j++){
                        dataArray=allData.get(j).split("\\t");
                        temp=dataArray[i];
                        //remove the "" in the csv file
                        //sometimes the length of the string is 0, for example, when there is no variable modification
                        if (temp.length()>0&&temp.charAt(0)=='"')
                            temp=temp.substring(1,temp.length()-1);
                        eachColumn.add(temp);
                    }
                    scanSummary.put(dataArrayIndex[i],eachColumn);
            }
        }
        
        catch(FileNotFoundException e){
            System.out.println("error1");
        }
        catch(IOException e){
       
            System.out.println("error2");
        }
    }
    public void findExpMSMS(File mgffile){
        expMSMS.clear();
        expMSMSBegin.clear();
        expMSMSEnd.clear();
        mgfFile=mgffile;
        scanOrder.clear();
        int nn;
        List<String> allDtaFile=scanSummary.get("Scan number");
         System.out.println(allDtaFile.get(0));
        try{
            BufferedReader MGFFile =new BufferedReader(new FileReader(mgfFile));
            String dataRow;
            //revised to get scan number instead of full name of dta file
            Pattern p=Pattern.compile("=.*\\.(\\d+)\\.\\d+\\.dta");
            Pattern p1=Pattern.compile("END IONS");
            Matcher m,m1;
            int lineNumber=0;
            String dtaFileName="";
            //System.out.println(dtaFileName);
            int n=0;
            String[] mzAndIntensity=new String[2];
            while ((dataRow=MGFFile.readLine())!=null){
                m=p.matcher(dataRow);
                if(m.find()){
                    n=n+1;
                    //System.out.println(dataRow);
                    //System.out.println(m.group(1));
                }
                m.reset();
                
                lineNumber=lineNumber+1;
                //only load the spectra listed in the scan summary file
                if(m.find()&&allDtaFile.indexOf(m.group(1))!=-1){
                    dtaFileName=m.group(1);
                    scanOrder.put(dtaFileName,n);
                    expMSMSBegin.put(dtaFileName,lineNumber+3);
                    MGFFile.readLine();
                    MGFFile.readLine();
                    lineNumber=lineNumber+2;
                    HashMap<Double,Double> oneSpectrum=new LinkedHashMap<>();
                    double preMz=0,curMz=0,preIn=0,curIn=0;
                    int indexes=0;
                    while((dataRow=MGFFile.readLine())!=null){
                        lineNumber=lineNumber+1;
                        m1=p1.matcher(dataRow);
                        if(!m1.find()){
                            mzAndIntensity=dataRow.split(" ");
                            curMz=Double.parseDouble(mzAndIntensity[0]);
                            curIn=Double.parseDouble(mzAndIntensity[1]);
                            /*if(curMz==preMz){
                                if(curIn>preIn){
                                    preIn=curIn;
                                    oneSpectrum.put(curMz,curIn);
                                }
                            }
                            else{
                                oneSpectrum.put(curMz,curIn);
                                preIn=curIn;
                            }
                            preMz=curMz;*/
                            if(oneSpectrum.get(curMz)==null){
                                oneSpectrum.put(curMz,curIn);
                                indexes=0;
                            }
                            else{
                                indexes=indexes+1;
                                oneSpectrum.put(curMz+0.00001*indexes,curIn);
                            }
                            
                        }
                        else{
                            expMSMSEnd.put(dtaFileName,lineNumber-1);
                            expMSMS.put(dtaFileName,oneSpectrum);
                            break;
                        }
                    }
                }
            }
            //System.out.println(expMSMS.size());
            //System.out.println(allDtaFile.size());
            
            if(expMSMS.size()<allDtaFile.size())
                System.out.println("scan number error");
            MGFFile.close();
            //System.out.println(expMSMSBegin);
            //System.out.println(expMSMS.get("bsa1_1hcdiw_10_01.5493.5493.3.dta"));
            //System.out.println(expMSMSEnd);
        }
        catch(FileNotFoundException e){
            ;
        }
        catch(IOException e){
            ;
        }
    }
    
}
