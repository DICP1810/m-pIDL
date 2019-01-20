/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package sixduplex;

import java.io.File;
import java.util.HashMap;
import java.util.LinkedHashMap;

/**
 *
 * @author shan
 */
public class RepeatedIon {
    double mzThreshold=400;
    double threshold=0.002;
    double setMassDiff=0.05844;
    public ReadScanSummary expPeptide;
    HashMap<String,HashMap<Double,Double>> ionPairs=new LinkedHashMap<>();
    HashMap<String,Double> ionRatios=new LinkedHashMap<>();
    public RepeatedIon(File scansummaryfile,File mgffile){
        expPeptide=new ReadScanSummary(scansummaryfile);
        expPeptide.findExpMSMS(mgffile);
    }
    public void calculateIonPair(){
       // for(String s:expMSMS.)
    
    }
}
