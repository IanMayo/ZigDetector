/**
 * Created by Bill on 1/23/2015.
 *A project to determine the Linear regression for maritime analytic using java
 * Modules such as apache commons maths libraries and Jfreechart are used for analysis and visualization
 */
import java.awt.FlowLayout;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.util.ArrayList;
import java.util.Date;
import java.util.Iterator;
import java.util.List;

import javax.swing.BoxLayout;
import javax.swing.JFrame;
import javax.swing.JPanel;

import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.MultiDirectionalSimplex;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.SimplexOptimizer;

public class ZigDetector {
	
	public static void main(String[] args) throws Exception {
		
		final String SCENARIO = "Scen1";
		
		Track ownshipTrack = new Track("data/" + SCENARIO +"_Ownship.csv");
		Track targetTrack = new Track("data/" + SCENARIO +"_Target.csv");
		Sensor sensor = new Sensor("data/" + SCENARIO +"_Sensor.csv");
		
		// create a holder for the data
		JFrame frame = new JFrame("Results");
		JPanel stack = new JPanel();
		stack.setLayout(new BoxLayout(stack, BoxLayout.Y_AXIS));
		frame.add(stack);
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

		
		// Now, we have to slice the data into ownship legs
		List<LegOfData> ownshipLegs = calculateLegs(ownshipTrack);
		
		// ok, time for the first plot
		Plotting.addOwnshipData(stack, "Ownship", ownshipTrack, ownshipLegs);
		Plotting.addOwnshipData(stack, "Target", targetTrack, null);
		
		frame.setSize(400, 800);
		frame.pack();
		
		long startTime = System.currentTimeMillis();
		
		// ok, work through the legs.  In the absence of a Discrete Optimisation algorithm we're taking a brue force approach.
		// Hopefully Craig can find an optimised alternative to this.
		for (Iterator<LegOfData> iterator = ownshipLegs.iterator(); iterator.hasNext();) {
			
			LegOfData thisLeg = (LegOfData) iterator.next();
			
			System.out.println(" handling leg:" + thisLeg);
			
			// ok, extract the relevant data
			List<Double> bearings = sensor.extractBearings(thisLeg.getStart(), thisLeg.getEnd());
			List<Long> times = sensor.extractTimes(thisLeg.getStart(), thisLeg.getEnd());
			
			// find the error score for the overall leg
	        MultivariateFunction wholeLeg = new ArcTanSolver(times, bearings); 
	        SimplexOptimizer wholeOptimizer = new SimplexOptimizer(1e-3, 1e-6); 
	        
	        int MAX_ITERATIONS = Integer.MAX_VALUE;
	        
			PointValuePair wholeLegOptimiser = wholeOptimizer.optimize( 
	                new MaxEval(MAX_ITERATIONS),
	                new ObjectiveFunction(wholeLeg), 
	                GoalType.MINIMIZE,
	                new InitialGuess(new double[] {bearings.get(0), 1, 1} ),//beforeBearings.get(0)
	                new MultiDirectionalSimplex(3)); 

			System.out.println(" whole leg score:" + wholeLegOptimiser.getValue().intValue());
			double[] key = wholeLegOptimiser.getKey();
			System.out.println("B:" + key[0] + " P:" + key[1] + " Q:" + key[2]);
			
			
			double bestScore = Double.MAX_VALUE;
			int bestIndex = -1;

			// make the two slices				
			final int BUFFER_REGION = 4; // the number of measurements to ignore whilst the target is turning 

			// how many points in this leg?
			int thisLegSize = times.size();
			int startIndex = 1 + BUFFER_REGION / 2;
			int endIndex = thisLegSize -1 - BUFFER_REGION / 2;
			
			for(int index=startIndex;index<endIndex;index++)
			{
				List<Long> theseTimes = times;
				List<Double> theseBearings = bearings;				
				
				// first the times
				List<Long> beforeTimes = theseTimes.subList(0, index - BUFFER_REGION / 2);
				List<Long> afterTimes = theseTimes.subList(index + BUFFER_REGION / 2, theseTimes.size()-1);
				
				// now the bearings
				List<Double> beforeBearings = theseBearings.subList(0, index);
				List<Double> afterBearings = theseBearings.subList(index, theseBearings.size()-1);
				
				
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
		        	        
		        // is this better?
		        if(sum < bestScore)
		        {
		        	// yes - store it.
		        	bestScore = sum;
		        	bestIndex = index;
		        }
			}			
			
	        System.out.println(" split sum:" + (int)bestScore + " at time " + new Date(times.get(bestIndex)));
	        
		}
		
		long elapsed = System.currentTimeMillis() - startTime;
		System.out.println("Elapsed:" + elapsed / 1000 + " secs");

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
		legs.add(new LegOfData("Ownship Leg 0"));
		
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
						legs.add(new LegOfData("Ownship Leg " + legs.size()));
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
