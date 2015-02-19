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
import java.awt.image.BufferedImage;
import java.beans.Transient;
import java.io.File;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.Iterator;
import java.util.List;

import javax.imageio.ImageIO;
import javax.swing.ImageIcon;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;

import junit.framework.TestCase;

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

	private static final int BUFFER_SIZE = 5;

	// how much RMS error we require on the Atan Curve before we
	// bother trying to slice the target leg
	private static final double ZIG_THRESHOLD = 0.005;

	// when to let the optimiser relax
	private static final double CONVERGE_TOLERANCE = 1e-6;

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

		LegStorer legStorer = new LegStorer();

		plotThis(inGrid, "Scen1", legStorer);
		// plotThis(inGrid, "Scen2a", legStorer);
		 plotThis(inGrid, "Scen2b", legStorer);
		  plotThis(inGrid, "Scen3", legStorer);
		 plotThis(inGrid, "Scen4", legStorer);

		if (inGrid.getComponentCount() == 1)
			grid.setColumns(1);

		frame.pack();

		long elapsed = System.currentTimeMillis() - startTime;
		System.out.println("Elapsed:" + elapsed / 1000 + " secs");

	}

	public static class LegStorer
	{
		private ArrayList<LegOfData> _legList;
		private TimeSeries _rmsScores;

		public void storeLeg(String scenario, long tStart, long tEnd, Sensor sensor, double rms)
		{
			System.out.println("Storing " + scenario + " : "
					+ dateF.format(new Date(tStart)) + " - "
					+ dateF.format(new Date(tEnd)));

			_legList.add(new LegOfData("Leg-" + (_legList.size()+1), tStart, tEnd));
			
			// create some RMS error scores
			List<Long> times = sensor.extractTimes(tStart, tEnd);
			for (Iterator<Long> iterator = times.iterator(); iterator.hasNext();)
			{
				Long long1 = (Long) iterator.next();
				_rmsScores.add(new FixedMillisecond(long1), rms);
			}
			
			
		}
		
		public void setRMSScores(TimeSeries series)
		{
			_rmsScores = series;
		}

		public void setLegList(ArrayList<LegOfData> legList)
		{
			_legList = legList;
		}

		public List<LegOfData> getLegs()
		{
			return _legList;
		}
	}

	public static void plotThis(Container container, String scenario,
			LegStorer legStorer) throws Exception
	{

		// load the data
		Track ownshipTrack = new Track("data/" + scenario + "_Ownship.csv");
		Track targetTrack = new Track("data/" + scenario + "_Target.csv");
		Sensor sensor = new Sensor("data/" + scenario + "_Sensor.csv");

		// slice the target legs, to help assess performance

		// Now, we have to slice the data into ownship legs
		// List<LegOfData> targetLegs = calculateLegs(targetTrack);
		List<LegOfData> ownshipLegs = identifyLegs(ownshipTrack);
		// ownshipLegs = ownshipLegs.subList(1, 2); // just play with the first leg

		// create the combined plot - where we show all our data
		CombinedDomainXYPlot combinedPlot = Plotting.createPlot();

		// get ready to store the results runs
		TimeSeriesCollection legResults = new TimeSeriesCollection();
		
		TimeSeries rmsScores = new TimeSeries("RMS Scores", FixedMillisecond.class);
		legStorer.setRMSScores(rmsScores);

		List<Long> turnMarkers = new ArrayList<Long>();
		legStorer.setLegList(new ArrayList<LegOfData>());
		
		// ok, work through the legs. In the absence of a Discrete Optimisation
		// algorithm we're taking a brue force approach.
		// Hopefully Craig can find an optimised alternative to this.
		for (Iterator<LegOfData> iterator = ownshipLegs.iterator(); iterator
				.hasNext();)
		{
			LegOfData thisLeg = (LegOfData) iterator.next();

			// ok, slice the data for this leg
			sliceThis(scenario, thisLeg.getStart(), thisLeg.getEnd(), sensor,
					legStorer);

			// create a placeholder for the overall score for this leg
			TimeSeries atanBar = new TimeSeries("ATan " + thisLeg.getName(),
					FixedMillisecond.class);
			legResults.addSeries(atanBar);

			// create a placeholder for the individual time slice experiments
			TimeSeries thisSeries = new TimeSeries(thisLeg.getName() + " Slices",
					FixedMillisecond.class);
			legResults.addSeries(thisSeries);

		}

		// show the track data (it contains the results)
		Plotting.plotSingleVesselData(combinedPlot, "O/S", ownshipTrack,
				ownshipLegs, new Color(0f, 0f, 1.0f, 0.2f), new Color(0f, 0f, 1.0f),
				null, timeEnd);

		Plotting.plotSingleVesselData(combinedPlot, "Tgt", targetTrack, legStorer
				.getLegs(), new Color(1.0f, 0f, 0f, 0.2f), new Color(1.0f, 0f, 0f),
				turnMarkers, timeEnd);

		Plotting.plotSensorData(combinedPlot, sensor.getTimes(), sensor.getBearings(), rmsScores);

		// insert the calculated P & Q
//		Plotting.plotPQData(combinedPlot, "Calculated", pqSeriesColl, null);

		// ok, also plot the leg attempts
//		Plotting.addLegResults(combinedPlot, legResults, turnMarkers);

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

	private static void sliceThis(String scenario, long curStart, long curEnd,
			Sensor sensor, LegStorer legStorer)
	{
		System.out.println("  Trying to slice : " + dateF.format(new Date(curStart))
				+ " - " + dateF.format(new Date(curEnd)) + " " + curStart + " to "
				+ curEnd);

		// ok, find the best slice
		// prepare the data
		List<Double> thisLegBearings = sensor.extractBearings(curStart, curEnd);
		List<Long> thisLegTimes = sensor.extractTimes(curStart, curEnd);
		
		if(thisLegBearings.size() == 0)
			return;

		Minimisation wholeLeg = optimiseThis(thisLegTimes, thisLegBearings,
				thisLegBearings.get(0));
		// is this slice acceptable?
		if (wholeLeg.getMinimum() < ZIG_THRESHOLD)
		{
			legStorer.storeLeg(scenario, curStart, curEnd, sensor, wholeLeg.getMinimum());
		}
		else
		{
			// ok, we'll have to slice it
			double bestScore = Double.MAX_VALUE;
			int bestSlice = -1;
			long sliceTime = -1;

			// find the optimal first slice
			for (int index = 0; index < thisLegTimes.size(); index++)
//			for (int index = 12; index < 17; index++)
			{
				// what's the total score for slicing at this index?
				double sum = sliceLeg(index, thisLegBearings, thisLegTimes, BUFFER_SIZE);

//				System.out.println("score for:" + dateF.format(new Date(thisLegTimes.get(index))) + " is " +  sum);
				
				// is this better?
				if ((sum > 0) && (sum < bestScore))
				{
					// yes - store it.
					bestScore = sum;
					bestSlice = index;
					sliceTime = thisLegTimes.get(index);
				}
			}
			
			// right, how did we get on?
			if (sliceTime != -1)
			{
				System.out.println("  Best slice at:" + dateF.format(new Date(sliceTime))
						+ " index:" + bestSlice
						+ " score:" + bestScore);

				// is this slice acceptable?
				if (bestScore < ZIG_THRESHOLD)
				{
					legStorer.storeLeg(scenario, curStart, sliceTime, sensor, bestScore);

					// have a look at the rest of the leg
					sliceThis(scenario, sliceTime + 60000, curEnd, sensor, legStorer);
				}
				else
				{
					List<Double> trimLegBearings = sensor.extractBearings(curStart, sliceTime);
					List<Long> trimLegTimes = sensor.extractTimes(curStart, sliceTime);
					
					// ok, see if we can reduce the buffer size
					boolean found = false;
					int bufferLen = 1;
					int maxBuffer = Math.min(trimLegTimes.size()-3, 7);
					
					
					while(!found && bufferLen < maxBuffer)
					{						
						List<Long> beforeTimes = trimLegTimes.subList(0, trimLegTimes.size() - bufferLen);
						List<Double> beforeBearings = trimLegBearings.subList(0, trimLegTimes.size() - bufferLen);
						
		//				System.out.println(" trimming:" + _outDates(beforeTimes));

						
						Minimisation beforeOptimiser = optimiseThis(beforeTimes, beforeBearings,
								beforeBearings.get(0));
						double sum = beforeOptimiser.getMinimum();
						
						if(sum < ZIG_THRESHOLD)
						{
							// ok, we can move on.
							found = true;
							legStorer.storeLeg(scenario, curStart, beforeTimes.get(beforeTimes.size()-1), sensor, sum);
							curStart = sliceTime + 1;
						}						

						bufferLen++;
					}
					
					// do we have enough left to bother with?
					sliceThis(scenario, curStart, curEnd, sensor, legStorer);
				}
			}
			else
			{
				System.out.println("slicing complete, can't slice");
			}
		}

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
	 * @param overallScore
	 *          the overall score for this leg
	 * @param BUFFER_REGION
	 * @param straightBar
	 * @param thisSeries
	 * @return
	 */
	private static double sliceLeg(int trialIndex, List<Double> bearings,
			List<Long> times, int bufferSize)
	{

		int legOneEnd = getEnd(0, times.size(), bufferSize, trialIndex);
		int legTwoStart = getStart(0, times.size(), bufferSize, trialIndex);

		List<Long> theseTimes = times;
		List<Double> theseBearings = bearings;

		Date thisD = new Date(times.get(trialIndex));

		// if((legOneEnd == -1) || (legTwoStart == -1))
		// return Double.MAX_VALUE;

		double beforeScore = 0;
		double afterScore = 0;

		@SuppressWarnings("unused")
		String msg = dateF.format(thisD);

		Minimisation beforeOptimiser = null;
		Minimisation afterOptimiser = null;

		if (legOneEnd != -1)
		{
			List<Long> beforeTimes = theseTimes.subList(0, legOneEnd);
			List<Double> beforeBearings = theseBearings.subList(0, legOneEnd);
			beforeOptimiser = optimiseThis(beforeTimes, beforeBearings,
					beforeBearings.get(0));
			beforeScore = beforeOptimiser.getMinimum();
			msg += " BEFORE:" + dateF.format(times.get(0)) + "-"
					+ dateF.format(times.get(legOneEnd)) + " ";
//			System.out.println(" before:" + _outDates(beforeTimes));
		}

		if (legTwoStart != -1)
		{
			List<Long> afterTimes = theseTimes.subList(legTwoStart,
					theseTimes.size() - 1);
			List<Double> afterBearings = theseBearings.subList(legTwoStart,
					theseTimes.size() - 1);
			afterOptimiser = optimiseThis(afterTimes, afterBearings,
					afterBearings.get(0));
			afterScore = afterOptimiser.getMinimum();
			msg += " AFTER:" + dateF.format(times.get(legTwoStart)) + "-"
					+ dateF.format(times.get(times.size() - 1)) + " ";
//			System.out.println(" after:" + _outDates(afterTimes));
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
	private static List<LegOfData> identifyLegs(Track track)
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
					if (legs.get(legs.size() - 1).initialised())
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
			double mean = runningSum / _times.size();

			double rms = Math.sqrt(mean);
			return rms;
		}

	}

	@SuppressWarnings("unused")
	private static double _splineErrorfor(List<Long> times, List<Double> bearings,
			Minimisation wholeLegOptimiser)
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
	
		// PolyTrendLine pt = new PolyTrendLine();
		// pt.setValues(y, x);
	
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
			// double thisForecast = pt.predict(thisMeasured);
			double thisError = Math.pow(thisForecast - thisMeasured, 2);
			runningSum += thisError;
		}
	
		for (int i = 0; i < times.size(); i++)
		{
			x[i] = (times.get(i) - startTime) / 1000d;
			y[i] = bearings.get(i);
			// System.out.println(x[i] + ", " + y[i] + "," + pt.predict(x[i]));
			System.out.println(x[i] + ", " + y[i] + "," + cs.interpolate(x[i]));
		}
	
		long start = times.get(0) / 1000;
		long lastTime = times.get(times.size() - 1) / 1000;
		for (long i = start; i < lastTime - 60; i += 60)
		{
			long thisT = i - start;
			// System.out.println(thisT + ",0 ," + pt.predict(thisT));
			System.out.println(thisT + ", ," + cs.interpolate(thisT));
		}
	
		return runningSum;
	}

	private static double calcForecast(double B, double P, double Q,
			double elapsedSecs)
	{
		double dX = Math.cos(Math.toRadians(B)) + Q * elapsedSecs;
		double dY = Math.sin(Math.toRadians(B)) + P * elapsedSecs;
		return Math.toDegrees(Math.atan2(dY, dX));
	}

	@SuppressWarnings("unused")
	private static TimeSeriesCollection _calculatePQ(Track ownship, Track target,
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

	@SuppressWarnings("unused")
	private static String _outDates(List<Long> times)
	{
		String res = dateF.format(times.get(0)) + "-"
				+ dateF.format(times.get(times.size() - 1));
		return res;
	}

	@SuppressWarnings("unused")
	private static double _splineErrorfor2(List<Long> times,
			List<Double> bearings, Minimisation wholeLegOptimiser)
	{
		long startTime = times.get(0);

		double[] keys = wholeLegOptimiser.getParamValues();
		double B = keys[0];
		double P = keys[1];
		double Q = keys[2];

		double runningSum = 0;
		for (int i = 0; i < times.size(); i++)
		{
			long thisT = (times.get(i) - startTime) / 1000;
			double thisB = bearings.get(i);
			double calcB = calcForecast(B, P, Q, thisT);

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
