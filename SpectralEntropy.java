
package moa.classifiers.core.driftdetection;


public class SpectralEntropy {

	private double[] powerSpectrumArray;
	
	private int slideWindowSize;
	
	//variables used into fft method.
    int swapPos;
    int bits;
    double[] swapBufferElement;
    double theta;
    double[] cusumPowerSpectrum;
    int evenIndex;
    int oddIndex;
    double[] even;
    double[] odd;
    double[] exp;
    
	
	public SpectralEntropy(int slideWindowSize) {
		
        // radix 2 Cooley-Tukey FFT
        if (slideWindowSize % 2 != 0) {
            throw new IllegalArgumentException("slideWindowSize is not a power of 2");
        }else {
        	this.slideWindowSize = slideWindowSize;
        }
		
        this.powerSpectrumArray = new double[this.slideWindowSize]; 
        
        this.even = new double[] {0.0,0.0};
        this.odd = new double[] {0.0,0.0};
        this.exp = new double[] {0.0,0.0};
        this.swapBufferElement = new double[] {0.0,0.0};
        this.cusumPowerSpectrum = new double[] {0.0,0.0};
        this.theta = 0.0;
        this.swapPos = 0;
        this.bits = 0;
        this.evenIndex = 0;
        this.oddIndex = 1;

        
	}
	
	public double run(double[][] timeSeries) {
		
        // radix 2 Cooley-Tukey FFT
        if (timeSeries.length  !=  this.slideWindowSize) {
            throw new IllegalArgumentException("Time series length is different of slide window size.");
        }
		        
        double cusumPowerSpectrum[] = {0.0,0.0};
        double spectralEntropy = 0.0;
        double auxNorm = 0.0;
 
        cusumPowerSpectrum = fft(timeSeries, this.powerSpectrumArray);
        
        //Compute the entropy using single-sided power spectrum (include DC).
        cusumPowerSpectrum[0] += cusumPowerSpectrum[1];
        auxNorm = this.powerSpectrumArray[0]/cusumPowerSpectrum[0];//DC
        spectralEntropy += auxNorm*Math.log(auxNorm);
        for (int i = 1; i <= powerSpectrumArray.length/2; i++) {
        	if (this.powerSpectrumArray[i] > 0.0) {
        		
        		//converting from a two-sided power spectrum to a single-sided power spectrum
        		auxNorm = 2*this.powerSpectrumArray[i]/cusumPowerSpectrum[0];
    			spectralEntropy += auxNorm*Math.log(auxNorm);
        	}
        }
        
        //Normalization
        spectralEntropy = -spectralEntropy*(1/Math.log(powerSpectrumArray.length/2));
		
		return spectralEntropy;
	}
	
	
	
	private int bitReverse(int n, int bits) {
        int reversedN = n;
        int count = bits - 1;
 
        n >>= 1;
        while (n > 0) {
            reversedN = (reversedN << 1) | (n & 1);
            count--;
            n >>= 1;
        }
 
        return ((reversedN << count) & ((1 << bits) - 1));
	 }
	 
 
	/**
	 * Compute the  Fast Fourier Transform of a time series and in addition return the power spectrum.
	 * Obtained and modified from: https://rosettacode.org/wiki/Fast_Fourier_transform
	 * 
	 * @param buffe:r a two dimensional array which the real part is the input sequence.
	 * @param powerSpectrum: a array to store the power spectrum energy.
	 */
	 private double[] fft(double[][] buffer, double[] powerSpectrum) {
		 
        // radix 2 Cooley-Tukey FFT
        if (buffer.length % 2 != 0) {
            throw new IllegalArgumentException("buffer size is not a power of 2.");
        }
        
        this.swapPos = 0;
        this.swapBufferElement[0] = 0.0;
        this.swapBufferElement[1] = 0.0;
        this.theta = 0.0;
        this.cusumPowerSpectrum[0] = 0.0;
        this.cusumPowerSpectrum[1] = 0.0;
        this.evenIndex = 0;
        this.oddIndex = 1;
        this.even[0] = 0.0;
        this.even[1] = 0.0;
        this.odd[0] = 0.0;
        this.odd[1] = 0.0;
        this.exp[0] = 0.0;
        this.exp[1] = 0.0;
		
        this.bits = (int) (Math.log(buffer.length) / Math.log(2));
        
        for (int j = 1; j < buffer.length / 2; j++) {
 
            swapPos = bitReverse(j, bits);
           
            swapBufferElement[0] = buffer[j][0];
            swapBufferElement[1] = buffer[j][1];

            buffer[j][0] = buffer[swapPos][0];
            buffer[j][1] = buffer[swapPos][1];

            buffer[swapPos][0] = swapBufferElement[0];
            buffer[swapPos][1] = swapBufferElement[1];
            
        }
 
        for (int N = 2; N <= buffer.length; N <<= 1) {
        	
            for (int i = 0; i < buffer.length; i += N) {
                for (int k = 0; k < N / 2; k++) {
 
                    evenIndex = i + k;
                    oddIndex = i + k + (N / 2);
                    
                    even[0] = buffer[evenIndex][0];
                    even[1] = buffer[evenIndex][1];
                    
                    odd[0] = buffer[oddIndex][0];
                    odd[1] = buffer[oddIndex][1];
                    
                    theta = (-2 * Math.PI * k) / (double) N;
                    exp[0] = Math.cos(theta) * odd[0] - Math.sin(theta) * odd[1];
                    exp[1] = Math.cos(theta) * odd[1] + Math.sin(theta) * odd[0];

                    buffer[evenIndex][0] = even[0] + exp[0];
                    buffer[evenIndex][1] = even[1] + exp[1];

                    buffer[oddIndex][0] = even[0] - exp[0];
                    buffer[oddIndex][1] = even[1] - exp[1];
                    
                    powerSpectrum[evenIndex] = Math.hypot(buffer[evenIndex][0], buffer[evenIndex][1]);
                    powerSpectrum[evenIndex] = powerSpectrum[evenIndex]*powerSpectrum[evenIndex];
                    
                    powerSpectrum[oddIndex] = Math.hypot(buffer[oddIndex][0], buffer[oddIndex][1]);
                    powerSpectrum[oddIndex] = powerSpectrum[oddIndex]*powerSpectrum[oddIndex];
                    
                    if ( N == buffer.length ) {
                    	cusumPowerSpectrum[0] += powerSpectrum[evenIndex];
                    	cusumPowerSpectrum[1] += powerSpectrum[oddIndex];
                    	
                    }
                    
                   
                }
            }
        }
        
        return cusumPowerSpectrum;
	 }
	 
	 

	 
	 
}
