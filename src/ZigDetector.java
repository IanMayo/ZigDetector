/**
 * Created by Bill on 1/23/2015.
 *A project to determine the Linear regression for maritime analytic using java
 * Modules such as apache commons maths libraries and Jfreechart are used for analysis and visualization
 */
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.GridLayout;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.beans.Transient;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.Iterator;
import java.util.List;

import javax.swing.JFrame;
import javax.swing.JPanel;

import junit.framework.TestCase;

import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.CombinedDomainXYPlot;
import org.jfree.data.time.FixedMillisecond;
import org.jfree.data.time.TimeSeries;
import org.jfree.data.time.TimeSeriesCollection;

import flanagan.interpolation.CubicSpline;
import flanagan.math.Minimisation;
import flanagan.math.MinimisationFunction;

public class ZigDetector
{

	private static final double CONVERGE_TOLERANCE = 1e-5;
	
//	final static long tgL1start = 1263297600000L; // 12:00:00 GMT
//	final static long osL1start = 1263297741000L; // 12:02:21 GMT
//	final static long tgL1end = 1263300091000L; // 12:41:31 GMT 2010
//	final static long tgL2start = 1263300279000L; // 12:44:39 GMT
//	final static long osL1end = 1263301172000L; // 12:59:32 GMT 2010
//	final static long osL2start = 1263301360000L; // 13:02:40 GMT
//	final static long tgL2end = 1263303616000L; // 13:40:16 GMT 2010
//	final static long tgL3start = 1263303804000L; // 13:43:24 GMT
//	final static long end = 1263304838000L; // 14:00:38 GMT 2010

	final static Long timeEnd = null; // osL1end;

	final static SimpleDateFormat dateF = new SimpleDateFormat("hh:mm:ss");
	final static DecimalFormat numF = new DecimalFormat(
			" 0000.0000000;-0000.0000000");

	public static void main(String[] args) throws Exception
	{

		// capture the start time (used for time elapsed at the end)
		long startTime = System.currentTimeMillis();

		// create a holder for the data
		final JFrame frame = createFrame();
		frame.setLocation(600, 50);
		Container container = frame.getContentPane();

		// ok, insert a grid
		JPanel inGrid = new JPanel();
		container.add(inGrid);
		GridLayout grid = new GridLayout(0, 2);
		inGrid.setLayout(grid);

		 plotThis(inGrid, "Scen1");
		plotThis(inGrid, "Scen2");
		 plotThis(inGrid, "Scen3");
		 plotThis(inGrid, "Scen4");

		if (inGrid.getComponentCount() == 1)
			grid.setColumns(1);

		frame.pack();

		long elapsed = System.currentTimeMillis() - startTime;
		System.out.println("Elapsed:" + elapsed / 1000 + " secs");

	}

	public static void plotThis(Container container, String scenario)
			throws Exception
	{

		// load the data
		Track ownshipTrack = new Track("data/" + scenario + "_Ownship.csv");
		Track targetTrack = new Track("data/" + scenario + "_Target.csv");
		Sensor sensor = new Sensor("data/" + scenario + "_Sensor.csv");

		// slice the target legs, to help assess performance
		List<LegOfData> targetLegs = calculateLegs(targetTrack);

		// Now, we have to slice the data into ownship legs
		List<LegOfData> ownshipLegs = calculateLegs(ownshipTrack);
//		ownshipLegs = ownshipLegs.subList(0, 1); // just play with the first leg

		// create the combined plot - where we show all our data
		CombinedDomainXYPlot combinedPlot = Plotting.createPlot();

		// ok create the plots of ownship & target tracks

		// get ready to store the results runs
		TimeSeriesCollection legResults = new TimeSeriesCollection();

		TimeSeriesCollection pqSeriesColl = calculatePQ(ownshipTrack, targetTrack,
				ownshipLegs);

		// get ready to store the fitted P & Q
		TimeSeries fittedBeforeP = new TimeSeries("P-bef-fit",
				FixedMillisecond.class);
		TimeSeries fittedBeforeQ = new TimeSeries("Q-bef-fit",
				FixedMillisecond.class);
		TimeSeries fittedAfterP = new TimeSeries("P-aft-fit",
				FixedMillisecond.class);
		TimeSeries fittedAfterQ = new TimeSeries("Q-aft-fit",
				FixedMillisecond.class);

		pqSeriesColl.addSeries(fittedBeforeP);
		pqSeriesColl.addSeries(fittedBeforeQ);
		pqSeriesColl.addSeries(fittedAfterP);
		pqSeriesColl.addSeries(fittedAfterQ);

		List<Long> valueMarkers = new ArrayList<Long>();

		// ok, work through the legs. In the absence of a Discrete Optimisation
		// algorithm we're taking a brue force approach.
		// Hopefully Craig can find an optimised alternative to this.
		for (Iterator<LegOfData> iterator = ownshipLegs.iterator(); iterator
				.hasNext();)
		{
			LegOfData thisLeg = (LegOfData) iterator.next();

			// ok, slice the data for this leg
			List<Double> bearings = sensor.extractBearings(thisLeg.getStart(),
					thisLeg.getEnd());
			List<Long> times = sensor.extractTimes(thisLeg.getStart(),
					thisLeg.getEnd());

			// find the error score for the overall leg
			Minimisation wholeLegOptimiser = optimiseThis(times, bearings,
					bearings.get(0));

			// find the RMS error for this arctan spline
			double splineError = splineErrorfor2(times, bearings, wholeLegOptimiser);

			System.out.println("&& Spline:" + splineError + " ArcTan:"
					+ wholeLegOptimiser.getMinimum());

			// look at the individual scores (though they're not of interest)
			// System.out.println("Whole Leg:" + thisLeg + " - " +
			// out(wholeLegOptimiser));

			// we will try to beat this score, so set as very high number
			double bestScore = Double.MAX_VALUE;
			int bestIndex = -1;

			Double overallScore = wholeLegOptimiser.getMinimum();

			// create a placeholder for the overall score for this leg
			TimeSeries atanBar = new TimeSeries("ATan " + thisLeg.getName(),
					FixedMillisecond.class);
			legResults.addSeries(atanBar);
			TimeSeries polyBar = new TimeSeries("RMS " + thisLeg.getName(),
					FixedMillisecond.class);
			legResults.addSeries(polyBar);
			
			// create a placeholder for the individual time slice experiments
			TimeSeries thisSeries = new TimeSeries(thisLeg.getName() + " Slices",
					FixedMillisecond.class);
			legResults.addSeries(thisSeries);

			// loop through the values in this leg
			// NOTE: this is a brute force algorithm - maybe we can find a
			// Discrete Optimisation equivalent
			for (int index = 0; index < times.size(); index++)
			{
				// what's the total score for slicing at this index?

				// if(index != 50)
				// continue;
				int legOneEnd = getEnd(0, times.size(), 5, index);
				int legTwoStart = getStart(0, times.size(), 5, index);

				// int legOneEnd = index - BUFFER_REGION / 2;
				// legOneEnd = Math.max(legOneEnd, 4);
				// int legTwoStart = legOneEnd + BUFFER_REGION;
				// legTwoStart = Math.min(legTwoStart, endIndex - 4);

				double sum = sliceLeg(legOneEnd, index, legTwoStart, bearings, times,
						fittedBeforeP, fittedBeforeQ, fittedAfterP, fittedAfterQ);

				thisSeries.add(new FixedMillisecond(times.get(index)), sum);
				atanBar.add(new FixedMillisecond(times.get(index)), overallScore);
				polyBar.add(new FixedMillisecond(times.get(index)), splineError);

				// is this better?
				if (sum < bestScore)
				{
					// yes - store it.
					bestScore = sum;
					bestIndex = index;
				}
			}

			// ok, decide if we should slice
			// System.out.println("leg: " + thisLeg.getName() + " whole:" +
			// overallScore + " best slice:" + bestScore);

			valueMarkers.add(times.get(bestIndex));

		}

		// show the track data (it contains the results)
		Plotting.plotSingleVesselData(combinedPlot, "O/S", ownshipTrack,
				ownshipLegs, new Color(0f, 0f, 1.0f, 0.2f), new Color(0f, 0f, 1.0f),
				null, timeEnd);

		Plotting.plotSingleVesselData(combinedPlot, "Tgt", targetTrack, targetLegs,
				new Color(1.0f, 0f, 0f, 0.2f), new Color(1.0f, 0f, 0f), valueMarkers,
				timeEnd);

		// insert the calculated P & Q
		// Plotting.plotPQData(combinedPlot, "Calculated", pqSeriesColl, null);

		// ok, also plot the leg attempts
		Plotting.addLegResults(combinedPlot, legResults, valueMarkers);


		// wrap the combined chart
		ChartPanel cp = new ChartPanel(new JFreeChart("Results for " + scenario
				+ " Tol:" + CONVERGE_TOLERANCE, JFreeChart.DEFAULT_TITLE_FONT,
				combinedPlot, true))
		{

			/**
					 * 
					 */
			private static final long serialVersionUID = 1L;

			@Override
			@Transient
			public Dimension getPreferredSize()
			{
				return new Dimension(700, 500);
			}

		};
		container.add(cp, BorderLayout.CENTER);

//		 BufferedImage wPic = ImageIO.read(new File("data/" + scenario +
//		 "_plot.png"));
//		 JLabel wIcon = new JLabel(new ImageIcon(wPic));
//		 container.add(wIcon, BorderLayout.SOUTH);

	}

	@SuppressWarnings("unused")
	private static double splineErrorfor(List<Long> times, List<Double> bearings, Minimisation wholeLegOptimiser)
	{
		long startTime = times.get(0);

		System.out.println("time, measured, spline");
		
		double[] x = new double[times.size()];
		double[] y = new double[times.size()];
		for (int i = 0; i < times.size(); i++)
		{
			x[i] = (times.get(i) - startTime) / 1000d;
			y[i] = bearings.get(i);
		}

		// ok, fit the spline
		CubicSpline cs = new CubicSpline(x, y);
		
	//	PolyTrendLine pt = new PolyTrendLine();
	//	pt.setValues(y, x);

		// and calculate the error sum
		double runningSum = 0;
		
		double[] keys = wholeLegOptimiser.getParamValues();

		double B = keys[0];
		double P = keys[1];
		double Q = keys[2];
		
		for (int i = 0; i < times.size(); i++)
		{
			long elapsedMillis = times.get(i) - startTime;
			double elapsedSecs = elapsedMillis / 1000d;
			
			double thisMeasured = bearings.get(i);
			double thisForecast = cs.interpolate(elapsedSecs);
//			double thisForecast = pt.predict(thisMeasured);
			double thisError = Math.pow(thisForecast - thisMeasured, 2);
			runningSum += thisError;
		}
		
		
		for (int i = 0; i < times.size(); i++)
		{
			x[i] = (times.get(i) - startTime) / 1000d;
			y[i] = bearings.get(i);
//			System.out.println(x[i] + ", " + y[i] + "," + pt.predict(x[i]));
			System.out.println(x[i] + ", " + y[i] + "," + cs.interpolate(x[i]));
		}
		
		long start = times.get(0)/ 1000;
		long lastTime = times.get(times.size()-1) / 1000;
		for(long i=start;i<lastTime-60;i+=60)
		{
			long thisT = i - start;
//			System.out.println(thisT + ",0 ," + pt.predict(thisT));
			System.out.println(thisT + ", ," + cs.interpolate(thisT));
		}
		

		return runningSum;
	}
	
	public interface TrendLine {
    public void setValues(double[] y, double[] x); // y ~ f(x)
    public double predict(double x); // get a predicted y for a given x
	}
	
	public static abstract class OLSTrendLine implements TrendLine {

		RealMatrix coef = null; // will hold prediction coefs once we get values

		protected abstract double[] xVector(double x); // create vector of values
														// from x

		protected abstract boolean logY(); // set true to predict log of y (note: y
											// must be positive)

		public void setValues(double[] y, double[] x) {
			if (x.length != y.length) {
				throw new IllegalArgumentException(String.format("The numbers of y and x values must be equal (%d != %d)",
						y.length, x.length));
			}
			double[][] xData = new double[x.length][];
			for (int i = 0; i < x.length; i++) {
				// the implementation determines how to produce a vector of
				// predictors from a single x
				xData[i] = xVector(x[i]);
			}
			if (logY()) { // in some models we are predicting ln y, so we replace
							// each y with ln y
				y = Arrays.copyOf(y, y.length); // user might not be finished with
												// the array we were given
				for (int i = 0; i < x.length; i++) {
					y[i] = Math.log(y[i]);
				}
			}
			OLSMultipleLinearRegression ols = new OLSMultipleLinearRegression();
			ols.setNoIntercept(true); // let the implementation include a constant
										// in xVector if desired
			ols.newSampleData(y, xData); // provide the data to the model
			coef = MatrixUtils.createColumnRealMatrix(ols.estimateRegressionParameters()); // get
																							// our
																							// coefs
		}

		public double predict(double x) {
			double yhat = coef.preMultiply(xVector(x))[0]; // apply coefs to xVector
			if (logY()) {
				yhat = (Math.exp(yhat)); // if we predicted ln y, we still need to
			}
			// get y
			return yhat;
		}
	}
	
	public static class PolyTrendLine extends OLSTrendLine {
	   @Override
	    protected double[] xVector(double x) {
	        return new double[]{1,x};
	    }

	    @Override
	    protected boolean logY() {return true;}
}
	

	private static TimeSeriesCollection calculatePQ(Track ownship, Track target,
			List<LegOfData> ownshipLegs)
	{
		TimeSeriesCollection ts = new TimeSeriesCollection();
		TimeSeries p = new TimeSeries("P-calc", FixedMillisecond.class);
		TimeSeries q = new TimeSeries("Q-calc", FixedMillisecond.class);
		ts.addSeries(p);
		ts.addSeries(q);

		// ok, now loop through
		long[] oTimes = ownship.getDates();
		double[] oCourses = ownship.getCourses();
		double[] oSpeeds = ownship.getSpeeds();
		double[] tCourses = target.getCourses();
		double[] tSpeeds = target.getSpeeds();

		double[] oX = ownship.getX();
		double[] oY = ownship.getY();
		double[] tX = target.getX();
		double[] tY = target.getY();

		Double R0 = null; // lMath.sqrt(Math.pow(oX[0], 2)+Math.pow(oY[0], 2));

		int thisLeg = 0;

		for (int i = 0; i < tY.length && thisLeg < ownshipLegs.size(); i++)
		{
			long thisT = oTimes[i];

			// ok, are we in a leg?
			if (thisT < ownshipLegs.get(thisLeg).getStart())
			{
				// ok, ignore
			}
			else
			{
				// see which leg we're in
				if (thisT < ownshipLegs.get(thisLeg).getEnd())
				{
					// hey, we're in the zone!
					// do we know range zero?
					if (R0 == null)
					{
						R0 = Math.sqrt(Math.pow(tX[i] - oX[i], 2)
								+ Math.pow(tY[i] - oY[i], 2));
					}

					// now for the other calcs
					double xO = oSpeeds[i] * Math.sin(Math.toRadians(oCourses[i]));
					double yO = oSpeeds[i] * Math.cos(Math.toRadians(oCourses[i]));
					double xT = tSpeeds[i] * Math.sin(Math.toRadians(tCourses[i]));
					double yT = tSpeeds[i] * Math.cos(Math.toRadians(tCourses[i]));

					double P = (xT - xO) / R0;
					double Q = (yT - yO) / R0;

					p.add(new FixedMillisecond(thisT), P);
					q.add(new FixedMillisecond(thisT), Q);
				}
				else
				{
					// ok, we've passed the end of the previous leg!
					thisLeg++;
					R0 = null;
				}
			}
		}

		return ts;
	}

	static Minimisation optimiseThis(List<Long> times, List<Double> bearings,
			double initialBearing)
	{
		// Create instance of Minimisation
		Minimisation min = new Minimisation();

		// Create instace of class holding function to be minimised
		FlanaganArctan funct = new FlanaganArctan(times, bearings);

		// initial estimates
		Double firstBearing = bearings.get(0);
		double[] start =
		{ firstBearing, 0.0D, 0.0D };

		// initial step sizes
		double[] step =
		{ 0.2D, 0.3D, 0.3D };

		// convergence tolerance
		double ftol = CONVERGE_TOLERANCE;

		// Nelder and Mead minimisation procedure
		min.nelderMead(funct, start, step, ftol);

		return min;
	}

	/**
	 * @param trialIndex
	 * @param bearings
	 * @param times
	 * @param fittedQ
	 * @param fittedP
	 * @param fittedAfterQ
	 * @param fittedAfterP
	 * @param overallScore
	 *          the overall score for this leg
	 * @param BUFFER_REGION
	 * @param straightBar
	 * @param thisSeries
	 * @return
	 */
	private static double sliceLeg(final int legOneEnd, int trialIndex,
			final int legTwoStart, List<Double> bearings, List<Long> times,
			TimeSeries fittedBeforeP, TimeSeries fittedBeforeQ,
			TimeSeries fittedAfterP, TimeSeries fittedAfterQ)
	{
		List<Long> theseTimes = times;
		List<Double> theseBearings = bearings;

		Date thisD = new Date(times.get(trialIndex));

		// if((legOneEnd == -1) || (legTwoStart == -1))
		// return Double.MAX_VALUE;

		double beforeScore = 0;
		double afterScore = 0;

		@SuppressWarnings("unused")
		String msg = dateF.format(thisD);

		if (legOneEnd != -1)
		{
			List<Long> beforeTimes = theseTimes.subList(0, legOneEnd);
			List<Double> beforeBearings = theseBearings.subList(0, legOneEnd);
			Minimisation beforeOptimiser = optimiseThis(beforeTimes, beforeBearings,
					beforeBearings.get(0));
			beforeScore = beforeOptimiser.getMinimum() / beforeTimes.size();
			msg += " BEFORE:" + dateF.format(times.get(0)) + "-"
					+ dateF.format(times.get(legOneEnd)) + " ";

			fittedBeforeP.add(new FixedMillisecond(times.get(trialIndex)),
					beforeOptimiser.getParamValues()[1]);
			fittedBeforeQ.add(new FixedMillisecond(times.get(trialIndex)),
					beforeOptimiser.getParamValues()[2]);
		}

		if (legTwoStart != -1)
		{
			List<Long> afterTimes = theseTimes.subList(legTwoStart,
					theseTimes.size() - 1);
			List<Double> afterBearings = theseBearings.subList(legTwoStart,
					theseTimes.size() - 1);
			Minimisation afterOptimiser = optimiseThis(afterTimes, afterBearings,
					afterBearings.get(0));
			afterScore = afterOptimiser.getMinimum() / afterTimes.size();
			msg += " AFTER:" + dateF.format(times.get(legTwoStart)) + "-"
					+ dateF.format(times.get(times.size() - 1)) + " ";

			fittedAfterP.add(new FixedMillisecond(times.get(trialIndex)),
					afterOptimiser.getParamValues()[1]);
			fittedAfterQ.add(new FixedMillisecond(times.get(trialIndex)),
					afterOptimiser.getParamValues()[2]);
		}

		// find the total error sum
		double sum = beforeScore + afterScore;

		// System.out.println(msg+ "SUM:" + sum);

		// DecimalFormat intF = new DecimalFormat("00");
		// System.out.println("index:"
		// + intF.format(trialIndex)
		// // + " time:" + times.get(trialIndex)
		// + " " + " Sum:" + numF.format(sum) + " index:"
		// + dateF.format(new Date(times.get(trialIndex))) + " before:"
		// + outDates(beforeTimes) + out(beforeOptimiser) + " num:"
		// + intF.format(beforeTimes.size()) + " after:" + outDates(afterTimes)
		// + out(afterOptimiser) + " num:" + intF.format(afterTimes.size()));

		return sum;
	}

	@SuppressWarnings("unused")
	private static String outDates(List<Long> times)
	{
		String res = dateF.format(times.get(0)) + "-"
				+ dateF.format(times.get(times.size() - 1));
		return res;
	}

	/**
	 * slice this data into ownship legs, where the course and speed are
	 * relatively steady
	 * 
	 * @param course_degs
	 * @param speed
	 * @param bearings
	 * @param elapsedTimes
	 * @return
	 */
	private static List<LegOfData> calculateLegs(Track track)
	{

		final double COURSE_TOLERANCE = 0.1; // degs / sec (just a guess!!)
		final double SPEED_TOLERANCE = 2; // knots / sec (just a guess!!)

		double lastCourse = 0;
		double lastSpeed = 0;
		long lastTime = 0;

		List<LegOfData> legs = new ArrayList<LegOfData>();
		legs.add(new LegOfData("Leg-1"));

		long[] times = track.getDates();
		double[] speeds = track.getSpeeds();
		double[] courses = track.getCourses();

		for (int i = 0; i < times.length; i++)
		{
			long thisTime = times[i];

			double thisSpeed = speeds[i];
			double thisCourse = courses[i];

			if (i > 0)
			{
				// ok, check out the course change rate
				double timeStepSecs = (thisTime - lastTime) / 1000;
				double courseRate = Math.abs(thisCourse - lastCourse) / timeStepSecs;
				double speedRate = Math.abs(thisSpeed - lastSpeed) / timeStepSecs;

				// are they out of range
				if ((courseRate < COURSE_TOLERANCE) && (speedRate < SPEED_TOLERANCE))
				{
					// ok, we're on a new leg - drop the current one
					legs.get(legs.size() - 1).add(thisTime);
				}
				else
				{
					// we may be in a turn. create a new leg, if we haven't done
					// so already
					if (legs.get(legs.size() - 1).size() != 0)
					{
						legs.add(new LegOfData("Leg-" + (legs.size() + 1)));
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

	private static class FlanaganArctan implements MinimisationFunction
	{
		final private List<Long> _times;
		final private List<Double> _bearings;

		public FlanaganArctan(List<Long> beforeTimes, List<Double> beforeBearings)
		{
			_times = beforeTimes;
			_bearings = beforeBearings;
		}

		// evaluation function
		public double function(double[] point)
		{
			double B = point[0];
			double P = point[1];
			double Q = point[2];

			double runningSum = 0;

			// ok, loop through the data
			for (int i = 0; i < _times.size(); i++)
			{
				long elapsedMillis = _times.get(i) - _times.get(0);
				double elapsedSecs = elapsedMillis / 1000d;
				double thisForecast = calcForecast(B, P, Q, elapsedSecs);
				double thisMeasured = _bearings.get(i);
				double thisError = Math.pow(thisForecast - thisMeasured, 2);
				runningSum += thisError;
			}
			return runningSum / _times.size();
		}

	}

	/**
	 * @return a frame to contain the results
	 */
	private static JFrame createFrame()
	{
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
		frame.setLayout(new BorderLayout());

		return frame;
	}

	public static String out(Minimisation res)
	{
		double[] key = res.getParamValues();
		String out = " B:" + numF.format(key[0]) + " P:" + numF.format(key[1])
				+ " Q:" + numF.format(key[2]) + " Sum:" + numF.format(res.getMinimum());

		return out;
	}

	public static int getEnd(int start, int end, int buffer, int index)
	{
		int res = -1;
		int MIN_SIZE = 3;
		int semiBuffer = buffer / 2;

		int earliestPossibleStartPoint = start + (MIN_SIZE - 1);

		if (index - semiBuffer > earliestPossibleStartPoint)
		{
			res = index - semiBuffer - 1;
		}
		return res;
	}

	public static int getStart(int start, int end, int buffer, int index)
	{
		int res = -1;
		int MIN_SIZE = 3;
		int semiBuffer = buffer / 2;

		int lastPossibleStartPoint = end - (MIN_SIZE - 1);

		if (index + semiBuffer < lastPossibleStartPoint)
		{
			res = index + semiBuffer + 1;
		}
		return res;
	}

	public static class TestMe extends TestCase
	{
		public void testStartTimes()
		{
			assertEquals(2, 5 / 2);

			assertEquals("correct", -1, getEnd(0, 15, 5, 0));
			assertEquals("correct", -1, getEnd(0, 15, 5, 1));
			assertEquals("correct", -1, getEnd(0, 15, 5, 2));
			assertEquals("correct", -1, getEnd(0, 15, 5, 3));
			assertEquals("correct", -1, getEnd(0, 15, 5, 4));
			assertEquals("correct", 2, getEnd(0, 15, 5, 5));
			assertEquals("correct", 3, getEnd(0, 15, 5, 6));
			assertEquals("correct", 4, getEnd(0, 15, 5, 7));
			assertEquals("correct", 5, getEnd(0, 15, 5, 8));
			assertEquals("correct", 6, getEnd(0, 15, 5, 9));
			assertEquals("correct", 7, getEnd(0, 15, 5, 10));
			assertEquals("correct", 8, getEnd(0, 15, 5, 11));
			assertEquals("correct", 9, getEnd(0, 15, 5, 12));
			assertEquals("correct", 10, getEnd(0, 15, 5, 13));
			assertEquals("correct", 11, getEnd(0, 15, 5, 14));
			assertEquals("correct", 12, getEnd(0, 15, 5, 15));

			assertEquals("correct", 3, getStart(0, 15, 5, 0));
			assertEquals("correct", 4, getStart(0, 15, 5, 1));
			assertEquals("correct", 5, getStart(0, 15, 5, 2));
			assertEquals("correct", 6, getStart(0, 15, 5, 3));
			assertEquals("correct", 7, getStart(0, 15, 5, 4));
			assertEquals("correct", 8, getStart(0, 15, 5, 5));
			assertEquals("correct", 9, getStart(0, 15, 5, 6));
			assertEquals("correct", 10, getStart(0, 15, 5, 7));
			assertEquals("correct", 11, getStart(0, 15, 5, 8));
			assertEquals("correct", 12, getStart(0, 15, 5, 9));
			assertEquals("correct", 13, getStart(0, 15, 5, 10));
			assertEquals("correct", -1, getStart(0, 15, 5, 11));
			assertEquals("correct", -1, getStart(0, 15, 5, 12));
			assertEquals("correct", -1, getStart(0, 15, 5, 13));
			assertEquals("correct", -1, getStart(0, 15, 5, 14));
			assertEquals("correct", -1, getStart(0, 15, 5, 15));

		}
	}

	private static double calcForecast(double B, double P, double Q, double elapsedSecs)
	{
		double dX = Math.cos(Math.toRadians(B)) + Q * elapsedSecs;
		double dY = Math.sin(Math.toRadians(B)) + P * elapsedSecs;
		return Math.toDegrees(Math.atan2(dY, dX));
	}

	private static double splineErrorfor2(List<Long> times, List<Double> bearings, Minimisation wholeLegOptimiser)
		{
			long startTime = times.get(0);
	
			double[] keys = wholeLegOptimiser.getParamValues();
			double B = keys[0];
			double P = keys[1];
			double Q = keys[2];
			
			double runningSum = 0;
			for(int i=0;i<times.size();i++)
			{
				long thisT = (times.get(i) - startTime) / 1000;
				double thisB = bearings.get(i);
				double calcB = calcForecast(B,  P,  Q, thisT);
				
				double error2 = Math.pow(calcB - thisB, 2);
				runningSum += error2;
			}
			
			// and the mean?
			double mean = runningSum / times.size();
			
			// lastly - sort out the root
			double root = Math.sqrt(mean);

			// done
			return root;
		}

}
