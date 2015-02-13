/**
 * Created by Bill on 1/23/2015.
 *A project to determine the Linear regression for maritime analytic using java
 * Modules such as apache commons maths libraries and Jfreechart are used for analysis and visualization
 */
import java.awt.BorderLayout;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.GridLayout;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.beans.Transient;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.Iterator;
import java.util.List;

import javax.swing.JFrame;
import javax.swing.JPanel;

import junit.framework.TestCase;

import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.CombinedDomainXYPlot;
import org.jfree.data.time.FixedMillisecond;
import org.jfree.data.time.TimeSeries;
import org.jfree.data.time.TimeSeriesCollection;

import flanagan.math.Minimisation;
import flanagan.math.MinimisationFunction;

public class ZigDetector
{

	final static long tgL1start = 1263297600000L; // 12:00:00 GMT
	final static long osL1start = 1263297741000L; // 12:02:21 GMT
	final static long tgL1end = 1263300091000L; // 12:41:31 GMT 2010
	final static long tgL2start = 1263300279000L; // 12:44:39 GMT
	final static long osL1end = 1263301172000L; // 12:59:32 GMT 2010
	final static long osL2start = 1263301360000L; // 13:02:40 GMT
	final static long tgL2end = 1263303616000L; // 13:40:16 GMT 2010
	final static long tgL3start = 1263303804000L; // 13:43:24 GMT
	final static long end = 1263304838000L; // 14:00:38 GMT 2010

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
		inGrid.setLayout(new GridLayout(1, 0));

		// plotThis(inGrid, "Scen1");
		plotThis(inGrid, "Scen2");

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

		// Now, we have to slice the data into ownship legs
		List<LegOfData> ownshipLegs = calculateLegs(ownshipTrack);
		// ownshipLegs = ownshipLegs.subList(1, 2); // just play with the first leg

		// create the combined plot - where we show all our data
		CombinedDomainXYPlot combinedPlot = Plotting.createPlot();

		// ok create the plots of ownship & target tracks
		Plotting.addOwnshipData(combinedPlot, "O/S ", ownshipTrack, ownshipLegs,
				null, timeEnd);

		// get ready to store the results runs
		TimeSeriesCollection legResults = new TimeSeriesCollection();

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

			// look at the individual scores (though they're not of interest)
			System.out.println("Whole Leg:" + out(wholeLegOptimiser));

			// we will try to beat this score, so set as very high number
			double bestScore = Double.MAX_VALUE;
			int bestIndex = -1;

			Double overallScore = wholeLegOptimiser.getMinimum();

			final int BUFFER_REGION = 2; // the number of measurements to ignore
			// whilst the target is turning

			// how many points in this leg?
			int thisLegSize = times.size();
			int startIndex = 2;
			int endIndex = thisLegSize - 3;

			// create a placeholder for the overall score for this leg
			TimeSeries straightBar = new TimeSeries("Whole " + thisLeg.getName(),
					FixedMillisecond.class);
			legResults.addSeries(straightBar);

			// create a placeholder for the individual time slice experiments
			TimeSeries thisSeries = new TimeSeries(thisLeg.getName(),
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

				double sum = sliceLeg(legOneEnd, index, legTwoStart, bearings, times);

				thisSeries.add(new FixedMillisecond(times.get(index)), sum);
				straightBar.add(new FixedMillisecond(times.get(index)), overallScore);

				// is this better?
				if (sum < bestScore)
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
		Plotting.addOwnshipData(combinedPlot, "Tgt ", targetTrack, null,
				valueMarkers, timeEnd);

		// wrap the combined chart
		ChartPanel cp = new ChartPanel(new JFreeChart("Results for " + scenario,
				JFreeChart.DEFAULT_TITLE_FONT, combinedPlot, true))
		{

			/**
					 * 
					 */
			private static final long serialVersionUID = 1L;

			@Override
			@Transient
			public Dimension getPreferredSize()
			{
				return new Dimension(1100, 800);
			}

		};
		container.add(cp, BorderLayout.CENTER);

		// System.exit(0);
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
		double[] start ={firstBearing, 0.0D, 0.0D};

		// initial step sizes
		double[] step =
		{ 0.2D, 0.6D, 0.2D };

		// convergence tolerance
		double ftol = 1e-7;

		// Nelder and Mead minimisation procedure
		min.nelderMead(funct, start, step, ftol);

		return min;
	}

	/**
	 * @param trialIndex
	 * @param bearings
	 * @param times
	 * @param overallScore
	 *          the overall score for this leg
	 * @param BUFFER_REGION
	 * @param straightBar
	 * @param thisSeries
	 * @return
	 */
	private static double sliceLeg(final int legOneEnd, int trialIndex,
			final int legTwoStart, List<Double> bearings, List<Long> times)
	{
		List<Long> theseTimes = times;
		List<Double> theseBearings = bearings;

		Date thisD = new Date(times.get(trialIndex));
		
//		if((legOneEnd == -1) || (legTwoStart == -1))
//				return Double.MAX_VALUE;
		
		double beforeScore = 0;
		double afterScore = 0;

		String msg = dateF.format(thisD);
		
		if (legOneEnd != -1)
		{
			List<Long> beforeTimes = theseTimes.subList(0, legOneEnd);
			List<Double> beforeBearings = theseBearings.subList(0, legOneEnd);
			Minimisation beforeOptimiser = optimiseThis(beforeTimes, beforeBearings,
					beforeBearings.get(0));
			beforeScore = beforeOptimiser.getMinimum() / beforeTimes.size();
			msg += " BEFORE:" + dateF.format(times.get(0))+"-"+dateF.format(times.get(legOneEnd)) + " ";
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
			msg += " AFTER:" + dateF.format(times.get(legTwoStart))+"-"+dateF.format(times.get(times.size()-1)) + " ";
		}

		// find the total error sum
		double sum = beforeScore + afterScore;
		
	//	System.out.println(msg+  "SUM:" + sum);

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

		private double calcForecast(double B, double P, double Q, double elapsedSecs)
		{
			double dX = Math.cos(Math.toRadians(B)) + Q * elapsedSecs;
			double dY = Math.sin(Math.toRadians(B)) + P * elapsedSecs;
			return Math.toDegrees(Math.atan2(dY, dX));
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

}
