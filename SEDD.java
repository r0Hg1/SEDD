
/**
 * Spectral Entropy Drift Detection Method.
 *
 * @author Roberto Souto Maior Barros (roberto@cin.ufpe.br)
 * @author Rohgi Chikushi (rtmc2@cin.ufpe.br)
 * @version $Revision: 1 $
 */

package moa.classifiers.core.driftdetection;

import moa.core.ObjectRepository;
import moa.options.FloatOption;
import moa.options.IntOption;
import moa.tasks.TaskMonitor;

public class SEDD extends AbstractChangeDetector {

	private static final long serialVersionUID = 1L;
	
    public IntOption numBitsBernoulliShiftExpansionOption = new IntOption(
            "numBitsBernoulliShiftExpansion",
            'i',
            "The number of bits for the finite expansion.",
            32, 4, 32);

    public IntOption minNumShiftOption = new IntOption(
            "minNumShift",
            'n',
            "The minimum number of shifts before permitting detecting change.",
            128, 32, 1024);
  
    public FloatOption magnitudeOption = new FloatOption(
    		"magnitude", 
    		'm',
            "Magnitude threshold", 
            0.37, 0.01, 1.0);
    
    private SpectralEntropy spectralEntropy;
   
    private double entropy;
    
    private double upperBound;
    private double lowerBound;
        
    private double pMagnitude;
    
    private long bernoulliShiftMask;
    private long bernoulliFactor;
    private long bernoulliShift;
    
    private int maxIndexData;
    private int indexBuffer;
    private int indexData;
    
    //private LimitedQueue<Double>  circularBufferShifts;
	private double[][] timeSeries;
	private double[] timeSeriesBuffer;
    
    private void initialize() {
    	
        this.upperBound = 0.0;
        this.lowerBound = 0.0;
        
        this.bernoulliFactor  = (long)Math.pow(2, this.numBitsBernoulliShiftExpansionOption.getValue());
    	this.bernoulliShiftMask = this.bernoulliFactor - 1;
        this.bernoulliShift = 0;
        
        this.indexData = 0;
        this.maxIndexData = this.minNumShiftOption.getValue() - 1;
    	
        this.pMagnitude = this.magnitudeOption.getValue();
    	
        this.timeSeries  = new double[this.maxIndexData + 1][2];
        this.timeSeriesBuffer = new double[this.maxIndexData + 1]; 
    	
    	this.spectralEntropy = new SpectralEntropy(this.maxIndexData + 1);
    	
    	this.resetLearning();
    }

    @Override
    public void resetLearning() {
    	this.bernoulliShift = 0;
    	this.indexData = 0;
    	
        this.upperBound = 0.0;
        this.lowerBound = 0.0;
            	
    	this.isChangeDetected = false;
    }

    @Override
    public void input(double predict) {
    	
        if (!this.isInitialized) {
            this.initialize();
            this.isInitialized = true;
        } else if (this.isChangeDetected) {
                resetLearning();  
        }
        
        this.bernoulliShift <<= 1;
        this.bernoulliShift += (long)predict;
        this.bernoulliShift &= this.bernoulliShiftMask;
        
      
    	//Insert a new value.
    	if (this.indexData > this.maxIndexData){
    		
    		this.entropy = this.spectralEntropy.run(this.timeSeries);
   		
        	if (this.entropy > 0.0 && this.entropy < 1.0) {
        		
            	if (this.upperBound < this.entropy ) {
        			this.upperBound = this.entropy;
        			
        			//Reset lowerBound
        			this.lowerBound = this.upperBound;
        		}
        		
        		if (this.lowerBound > this.entropy  ) {
        			
        			this.lowerBound = this.entropy;
        		}	
        		
        	}
    		
    	    if ( (this.upperBound - this.lowerBound) >= this.pMagnitude)  {
    	    	
    	    	this.isChangeDetected = true;
    	    	
    	    }
    	    
    	    for (indexBuffer = 1; indexBuffer <= this.maxIndexData; indexBuffer++) {
    	    	this.timeSeriesBuffer[indexBuffer - 1] = this.timeSeriesBuffer[indexBuffer];
    	    	this.timeSeries[indexBuffer - 1][0] = this.timeSeriesBuffer[indexBuffer];
    	    	this.timeSeries[indexBuffer - 1][1] = 0.0;
    	    	
    	    }
    	    
    	    this.indexData = this.maxIndexData;
    	    
        }
    	
        this.timeSeries[this.indexData][0] = (this.bernoulliShift/(double)this.bernoulliFactor);
        this.timeSeries[this.indexData][1] = 0.0;
        this.timeSeriesBuffer[this.indexData] = this.timeSeries[this.indexData][0]; 
        
        this.indexData++;
		
    }

    @Override
    public void getDescription(StringBuilder sb, int indent) {
        // TODO Auto-generated method stub
    }

    @Override
    protected void prepareForUseImpl(TaskMonitor monitor,
            ObjectRepository repository) {
        // TODO Auto-generated method stub
    }
}
