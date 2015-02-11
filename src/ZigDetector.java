/**
 * Created by Bill on 1/23/2015.
 *A project to determine the Linear regression for maritime analytic using java
 * Modules such as apache commons maths libraries and Jfreechart are used for analysis and visualization
 */
import java.awt.Dimension;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.beans.Transient;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import javax.swing.JFrame;

import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.MultiDirectionalSimplex;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.SimplexOptimizer;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.CombinedDomainXYPlot;
import org.jfree.data.time.FixedMillisecond;
import org.jfree.data.time.TimeSeries;
import org.jfree.data.time.TimeSeriesCollection;

public class ZigDetector {
	
	public static void main(String[] args) throws Exception {
		
		// declare which scenario we're running
		final String SCENARIO = "Scen1";
		
		// load the data
		Track ownshipTrack = new Track("data/" + SCENARIO +"_Ownship.csv");
		Track targetTrack = new Track("data/" + SCENARIO +"_Target.csv");
		Sensor sensor = new Sensor("data/" + SCENARIO +"_Sensor.csv");
		
		// create a holder for the data
		final JFrame frame = createFrame();

		// Now, we have to slice the data into ownship legs
		List<LegOfData> ownshipLegs = calculateLegs(ownshipTrack);
		
		// create the combined plot - where we show all our data
		CombinedDomainXYPlot combinedPlot = Plotting.createPlot();
		
		// ok create the plots of ownship & target tracks
		Plotting.addOwnshipData(combinedPlot, "O/S", ownshipTrack, ownshipLegs, null);
		
		// capture the start time (used for time elapsed at the end)
		long startTime = System.currentTimeMillis();
		
		// get ready to store the results runs
		TimeSeriesCollection legResults = new TimeSeriesCollection();
		
		List<Long> valueMarkers = new ArrayList<Long>();
		
		// ok, work through the legs.  In the absence of a Discrete Optimisation algorithm we're taking a brue force approach.
		// Hopefully Craig can find an optimised alternative to this.
		for (Iterator<LegOfData> iterator = ownshipLegs.iterator(); iterator.hasNext();) {
			
			LegOfData thisLeg = (LegOfData) iterator.next();
			
			// ok, slice the data for this leg
			List<Double> bearings = sensor.extractBearings(thisLeg.getStart(), thisLeg.getEnd());
			List<Long> times = sensor.extractTimes(thisLeg.getStart(), thisLeg.getEnd());
			
			// find the error score for the overall leg
	        MultivariateFunction wholeLeg = new ArcTanSolver(times, bearings);
	        
	        // 
	        SimplexOptimizer wholeOptimizer = new SimplexOptimizer(1e-3, 1e-6); 
	        
	        //
	        int MAX_ITERATIONS = Integer.MAX_VALUE;
	        
	        // calculate the overall score for this leg
			PointValuePair wholeLegOptimiser = wholeOptimizer.optimize( 
	                new MaxEval(MAX_ITERATIONS),
	                new ObjectiveFunction(wholeLeg), 
	                GoalType.MINIMIZE,
	                new InitialGuess(new double[] {bearings.get(0), 1, 1} ),//beforeBearings.get(0)
	                new MultiDirectionalSimplex(3)); 

			// look at the individual scores (though they're not of interest)
//			double[] key = wholeLegOptimiser.getKey();
//			System.out.println("B:" + key[0] + " P:" + key[1] + " Q:" + key[2]);
	
			// we will try to beat this score, so set as very high number
			double bestScore = Double.MAX_VALUE;
			int bestIndex = -1;

			final int BUFFER_REGION = 4; // the number of measurements to ignore whilst the target is turning 

			// how many points in this leg?
			int thisLegSize = times.size();
			int startIndex = 1 + BUFFER_REGION / 2;
			int endIndex = thisLegSize -1 - BUFFER_REGION / 2;

			// create a placeholder for the overall score for this leg
			TimeSeries straightBar = new TimeSeries("Whole " + thisLeg.getName(), FixedMillisecond.class);
			legResults.addSeries(straightBar);
			
			// create a placeholder for the individual time slice experiments
			TimeSeries thisSeries = new TimeSeries(thisLeg.getName(), FixedMillisecond.class);
			legResults.addSeries(thisSeries);
						
			// loop through the values in this leg
			// NOTE: this is a brute force algorithm - maybe we can find a Discrete Optimisation equivalent
			for(int index=startIndex;index<endIndex;index++)
			{
				// what's the total score for slicing at this index?
				double sum = sliceLeg(index, bearings, times, MAX_ITERATIONS,
						wholeLegOptimiser.getValue(), BUFFER_REGION, straightBar,
						thisSeries);
		        	        
		        // is this better?
		        if(sum < bestScore)
		        {
		        	// yes - store it.
		        	bestScore = sum;
		        	bestIndex = index;
		        }
			}			
			
			valueMarkers.add(times.get(bestIndex));
			
		}
		
		// ok, also plot the leg attempts
		Plotting.addLegResults(combinedPlot, legResults, valueMarkers);
		
		// show the target track (it contains the results)
		Plotting.addOwnshipData(combinedPlot, "Tgt ", targetTrack, null, valueMarkers);

		// wrap the combined chart 
		ChartPanel cp = new ChartPanel(new JFreeChart("Results for " + SCENARIO,
                JFreeChart.DEFAULT_TITLE_FONT, combinedPlot, true)){

					@Override
					@Transient
					public Dimension getPreferredSize() {
						return new Dimension(1000,800);
					}
			
		};
		frame.add(cp);		
		frame.pack();
		
		long elapsed = System.currentTimeMillis() - startTime;
		System.out.println("Elapsed:" + elapsed / 1000 + " secs");
	}

	/**
	 * @return a frame to contain the results
	 */
	private static JFrame createFrame() {
		JFrame frame = new JFrame("Results");
		frame.pack();
		frame.setVisible(true);
		frame.addWindowListener(new WindowAdapter()
        {
            @Override
            public void windowClosing(WindowEvent e)
            {
                System.out.println("Closed");
                e.getWindow().dispose();
            }
        });

		return frame;
	}

	/**
	 * @param trialIndex
	 * @param bearings
	 * @param times
	 * @param MAX_ITERATIONS
	 * @param overallScore the overall score for this leg
	 * @param BUFFER_REGION
	 * @param straightBar
	 * @param thisSeries
	 * @return
	 */
	private static double sliceLeg(int trialIndex, List<Double> bearings,
			List<Long> times, int MAX_ITERATIONS,
			double overallScore, final int BUFFER_REGION,
			TimeSeries straightBar, TimeSeries thisSeries) {
		List<Long> theseTimes = times;
		List<Double> theseBearings = bearings;				
		
		// first the times
		List<Long> beforeTimes = theseTimes.subList(0, trialIndex - BUFFER_REGION / 2);
		List<Long> afterTimes = theseTimes.subList(trialIndex + BUFFER_REGION / 2, theseTimes.size()-1);
		
		// now the bearings
		List<Double> beforeBearings = theseBearings.subList(0, trialIndex);
		List<Double> afterBearings = theseBearings.subList(trialIndex, theseBearings.size()-1);
		
		MultivariateFunction beforeF = new ArcTanSolver(beforeTimes, beforeBearings); 
		MultivariateFunction afterF = new ArcTanSolver(afterTimes, afterBearings); 
		    		
		SimplexOptimizer optimizerMult = new SimplexOptimizer(1e-3, 1e-6); 
		
		// fit a curve to the period before the turn
		PointValuePair beforeOptimiser = optimizerMult.optimize( 
		        new MaxEval(MAX_ITERATIONS),
		        new ObjectiveFunction(beforeF), 
		        GoalType.MINIMIZE,
		        new InitialGuess(new double[] {beforeBearings.get(0), 1, 1} ),//beforeBearings.get(0)
		        new MultiDirectionalSimplex(3)); 
		
		// fit a curve to the period after the turn
		PointValuePair afterOptimiser = optimizerMult.optimize( 
		        new MaxEval(MAX_ITERATIONS),
		        new ObjectiveFunction(afterF), 
		        GoalType.MINIMIZE,
		        new InitialGuess(new double[] {afterBearings.get(0), 1, 1} ),//afterBearings.get(0)
		        new MultiDirectionalSimplex(3)); 

		// find the total error sum
		double sum = beforeOptimiser.getValue() + afterOptimiser.getValue();
		
		thisSeries.add(new FixedMillisecond(times.get(trialIndex)), sum);
		straightBar.add(new FixedMillisecond(times.get(trialIndex)),  overallScore);
		return sum;
	}

	/** slice this data into ownship legs, where the course and speed are relatively steady
	 * 
	 * @param course_degs
	 * @param speed
	 * @param bearings
	 * @param elapsedTimes
	 * @return
	 */
	private static List<LegOfData> calculateLegs(Track track) {
		
		final double COURSE_TOLERANCE = 0.1;   // degs / sec (just a guess!!)
		final double SPEED_TOLERANCE = 2;   // knots / sec  (just a guess!!)
		
		double lastCourse = 0;
		double lastSpeed = 0;
		long lastTime = 0;
		
		List<LegOfData> legs = new ArrayList<LegOfData>();
		legs.add(new LegOfData("Leg-1"));
		
		long[] times = track.getDates();
		double[] speeds = track.getSpeeds();
		double[] courses = track.getCourses();
		
		for (int i = 0; i < times.length; i++) {
			long thisTime = times[i];
			
			double thisSpeed = speeds[i];
			double thisCourse = courses[i];
			
			if(i > 0)
			{
				// ok, check out the course change rate
				double timeStepSecs = (thisTime - lastTime)/1000;
				double courseRate = Math.abs(thisCourse - lastCourse) / timeStepSecs; 
				double speedRate = Math.abs(thisSpeed - lastSpeed) / timeStepSecs;
				
				// are they out of range
				if((courseRate < COURSE_TOLERANCE) && (speedRate < SPEED_TOLERANCE))
				{
					// ok, we're on a new leg - drop the current one
					legs.get(legs.size()-1).add(thisTime);
				}
				else
				{
					// we may be in a turn. create a new leg, if we haven't done so already
					if(legs.get(legs.size()-1).size() != 0)
					{
						legs.add(new LegOfData("Leg-" + (legs.size()+1)));
					}
				}
			}
			
			// ok, store the values
			lastTime = thisTime;
			lastCourse = thisCourse;
			lastSpeed = thisSpeed;
			
		}
		
		return legs;
	}

	/** function to generate sum of squares of errors for a single permutation of B,P,Q
	 * 
	 */
	private static class ArcTanSolver implements MultivariateFunction
	{
		final private List<Long> _times;
		final private List<Double> _bearings;

		public ArcTanSolver(List<Long> beforeTimes, List<Double> beforeBearings)
		{
			_times = beforeTimes;
			_bearings = beforeBearings;
		}

		@Override
		public double value(double[] point) {
        	double B = point[0];
        	double P = point[1];
        	double Q = point[2];
        	
        	double runningSum = 0;
        	
        	// ok, loop through the data
        	for (int i = 0; i < _times.size(); i++) {
        		Long elapsedSecs = _times.get(i);
				double thisForecast = calcForecast(B, P, Q, elapsedSecs / 1000d); 
				Double thisMeasured = _bearings.get(i);
				double thisError = Math.pow(thisForecast - thisMeasured, 2);
				runningSum += thisError;
			}        	

        //	System.out.println("B:" + (int)B + " P:" + (int)P + " Q:" + (int)Q + " sum:" + runningSum);
        	
            return runningSum;
		}

		private double calcForecast(double B, double P, double Q,
				Double elapsedSecs) {
//			return Math.toDegrees(Math.atan2(Math.sin(Math.toRadians(B))+P*elapsedSecs,Math.cos(Math.toRadians(B))+Q*elapsedSecs));
			return Math.toDegrees(Math.atan2(Math.cos(Math.toRadians(B))+Q*elapsedSecs, Math.sin(Math.toRadians(B))+P*elapsedSecs));
		}		
	}
}
