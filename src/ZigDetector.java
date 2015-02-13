/**
 * Created by Bill on 1/23/2015.
 *A project to determine the Linear regression for maritime analytic using java
 * Modules such as apache commons maths libraries and Jfreechart are used for analysis and visualization
 */
import java.awt.Dimension;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.beans.Transient;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.Iterator;
import java.util.List;

import javax.swing.JFrame;

import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.optim.ConvergenceChecker;
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

	final static long tgL1start = 1263297600000L; // 12:00:00 GMT
	final static long osL1start = 1263297741000L; // 12:02:21 GMT
	final static long tgL1end = 1263300091000L; // 12:41:31 GMT 2010
	final static long tgL2start = 1263300279000L; // 12:44:39 GMT
	final static long osL1end = 1263301172000L; // 12:59:32 GMT 2010
	final static long osL2start = 1263301360000L; // 13:02:40 GMT
	final static long tgL2end = 1263303616000L; // 13:40:16 GMT 2010
	final static long tgL3start = 1263303804000L; // 13:43:24 GMT
	final static long end = 1263304838000L; // 14:00:38 GMT 2010

	final static Long timeEnd =  osL1end;
	
	final static int MAX_ITERATIONS = Integer.MAX_VALUE;

	final static SimpleDateFormat dateF = new SimpleDateFormat("hh:mm:ss");
	final static DecimalFormat numF = new DecimalFormat(" 0000.0000000;-0000.0000000");

	final static int doTest = 0;

	public static void main2(String[] args) throws IOException {
		// declare which scenario we're running
		final String SCENARIO = "Scen1";

		// load the data
		// Track ownshipTrack = new Track("data/" + SCENARIO + "_Ownship.csv");
		// Track targetTrack = new Track("data/" + SCENARIO +"_Target.csv");
		Sensor sensor = new Sensor("data/" + SCENARIO + "_Sensor.csv");

		// ok, slice the data for this leg
		debugOptimiseThis("leg one", sensor, osL1start, osL1end);

		System.out.println("slice at:" + new Date(tgL1end));
		// debugOptimiseThis("leg one befor", sensor, osL1start, tgL1end);
		// debugOptimiseThis("leg one after", sensor, tgL2start, osL1end);

	//	debugOptimiseThis("test run", sensor, 1263300185000L, 1263301078000L);

		debugOptimiseThis("36 ", sensor, 1263297788000L, 1263299386000L);
		debugOptimiseThis("37 ", sensor, 1263297788000L, 1263299433000L);
		debugOptimiseThis("38 ", sensor, 1263297788000L, 1263299480000L);
		debugOptimiseThis("39 ", sensor, 1263297788000L, 1263299527000L);
		debugOptimiseThis("40 ", sensor, 1263297788000L, 1263299574000L);
		debugOptimiseThis("41 ", sensor, 1263297788000L, 1263299621000L);
		debugOptimiseThis("42 ", sensor, 1263297788000L, 1263299668000L);
		debugOptimiseThis("43 ", sensor, 1263297788000L, 1263299715000L);
		debugOptimiseThis("44 ", sensor, 1263297788000L, 1263299762000L);

	
		
//		index:36  Sum:0024.69 index:12:31:20 before:12:03:08-12:29:46 B:-0152.85 P:-0000.00 Q:-0000.00 Sum: 0000.21 after:12:32:07-12:57:58 B:-0152.70 P: 0000.00 Q: 0000.00 Sum: 0024.48
//		index:37  Sum:0006.54 index:12:32:07 before:12:03:08-12:30:33 B:-0152.85 P:-0000.00 Q:-0000.00 Sum: 0000.23 after:12:32:54-12:57:58 B:-0151.65 P: 0000.00 Q: 0000.00 Sum: 0006.31
//		index:38  Sum:0034.55 index:12:32:54 before:12:03:08-12:31:20 B:-0152.97 P:-0000.00 Q:-0000.00 Sum: 0000.02 after:12:33:41-12:57:58 B:-0153.25 P: 0000.00 Q: 0000.00 Sum: 0034.53
//		index:39  Sum:0055.86 index:12:33:41 before:12:03:08-12:32:07 B:-0153.01 P:-0000.00 Q:-0000.00 Sum: 0000.01 after:12:34:28-12:57:58 B:-0153.86 P: 0000.00 Q: 0000.00 Sum: 0055.85
//		index:40  Sum:0031.91 index:12:34:28 before:12:03:08-12:32:54 B:-0152.78 P:-0000.00 Q:-0000.00 Sum: 0000.56 after:12:35:15-12:57:58 B:-0153.30 P: 0000.00 Q: 0000.00 Sum: 0031.35
//		index:41  Sum:0002.35 index:12:35:15 before:12:03:08-12:33:41 B:-0152.77 P:-0000.00 Q:-0000.00 Sum: 0000.62 after:12:36:02-12:57:58 B:-0151.84 P: 0000.00 Q: 0000.00 Sum: 0001.74
//		index:42  Sum:0050.53 index:12:36:02 before:12:03:08-12:34:28 B:-0153.01 P:-0000.00 Q:-0000.00 Sum: 0000.01 after:12:36:49-12:57:58 B:-0154.01 P: 0000.00 Q: 0000.00 Sum: 0050.53
//		index:43  Sum:0001.26 index:12:36:49 before:12:03:08-12:35:15 B:-0153.01 P:-0000.00 Q:-0000.00 Sum: 0000.01 after:12:37:36-12:57:58 B:-0151.92 P: 0000.00 Q: 0000.00 Sum: 0001.25
//		index:44  Sum:0000.77 index:12:37:36 before:12:03:08-12:36:02 B:-0153.01 P:-0000.00 Q:-0000.00 Sum: 0000.01 after:12:38:23-12:57:58 B:-0151.91 P: 0000.00 Q: 0000.00 Sum: 0000.76
	}

	/**
	 * @param sensor
	 * @param start
	 * @param end
	 * @return
	 */
	private static PointValuePair debugOptimiseThis(String name, Sensor sensor,
			long start, long end) {
		List<Double> bearings = sensor.extractBearings(start, end);
		List<Long> times = sensor.extractTimes(start, end);

		PointValuePair wholeLegOptimiser = optimiseThis(times, bearings, bearings.get(0));

		// look at the individual scores (though they're not of interest)
		double[] key = wholeLegOptimiser.getKey();
		System.out.println(name + " times:" + outDates(times) + " B:"
				+ numF.format(key[0]) + " P:" + numF.format(key[1]) + " Q:"
				+ numF.format(key[2]) + " Sum:" + wholeLegOptimiser.getValue());

		return wholeLegOptimiser;
	}

	@SuppressWarnings("unused")
	public static void main(String[] args) throws Exception {

		if (doTest > 0) {
			main2(args);
			 System.exit(0);
		}

		System.out.println("============");

		// declare which scenario we're running
		final String SCENARIO = "Scen1";

		// load the data
		Track ownshipTrack = new Track("data/" + SCENARIO + "_Ownship.csv");
		Track targetTrack = new Track("data/" + SCENARIO + "_Target.csv");
		Sensor sensor = new Sensor("data/" + SCENARIO + "_Sensor.csv");

		// create a holder for the data
		final JFrame frame = createFrame();

		// Now, we have to slice the data into ownship legs
		List<LegOfData> ownshipLegs = calculateLegs(ownshipTrack);
		
		// just play with the first leg
		ownshipLegs = ownshipLegs.subList(0, 1);

		// create the combined plot - where we show all our data
		CombinedDomainXYPlot combinedPlot = Plotting.createPlot();

		// ok create the plots of ownship & target tracks
		Plotting.addOwnshipData(combinedPlot, "O/S ", ownshipTrack,
				ownshipLegs, null, timeEnd);

		// capture the start time (used for time elapsed at the end)
		long startTime = System.currentTimeMillis();

		// get ready to store the results runs
		TimeSeriesCollection legResults = new TimeSeriesCollection();

		List<Long> valueMarkers = new ArrayList<Long>();

		// ok, work through the legs. In the absence of a Discrete Optimisation
		// algorithm we're taking a brue force approach.
		// Hopefully Craig can find an optimised alternative to this.
		 for (Iterator<LegOfData> iterator = ownshipLegs.iterator(); iterator
		 .hasNext();) {
		 LegOfData thisLeg = (LegOfData) iterator.next();

	//	{
	//		LegOfData thisLeg = ownshipLegs.get(0);

			// ok, slice the data for this leg
			List<Double> bearings = sensor.extractBearings(thisLeg.getStart(),
					thisLeg.getEnd());
			List<Long> times = sensor.extractTimes(thisLeg.getStart(),
					thisLeg.getEnd());

			// find the error score for the overall leg
			PointValuePair wholeLegOptimiser = optimiseThis(times, bearings, bearings.get(0));

			// look at the individual scores (though they're not of interest)
			System.out.println("Whole Leg:" + out(wholeLegOptimiser));
			// double[] key = wholeLegOptimiser.getKey();
			// System.out.println(thisLeg + " B:" + (int) key[0] + " P:" +
			// key[1]
			// + " Q:" + key[2]);

			// we will try to beat this score, so set as very high number
			double bestScore = Double.MAX_VALUE;
			int bestIndex = -1;

			Double overallScore = wholeLegOptimiser.getValue();

			final int BUFFER_REGION = 2; // the number of measurements to ignore
											// whilst the target is turning

			// how many points in this leg?
			int thisLegSize = times.size();
			int startIndex = 2;
			int endIndex = thisLegSize - 3;

			// create a placeholder for the overall score for this leg
			TimeSeries straightBar = new TimeSeries("Whole "
					+ thisLeg.getName(), FixedMillisecond.class);
			legResults.addSeries(straightBar);

			// create a placeholder for the individual time slice experiments
			TimeSeries thisSeries = new TimeSeries(thisLeg.getName(),
					FixedMillisecond.class);
			legResults.addSeries(thisSeries);

			// loop through the values in this leg
			// NOTE: this is a brute force algorithm - maybe we can find a
			// Discrete Optimisation equivalent
			for (int index = startIndex; index < endIndex; index++) {
				// what's the total score for slicing at this index?

				// if(index != 50)
				// continue;

				int legOneEnd = index - BUFFER_REGION / 2;
				legOneEnd = Math.max(legOneEnd, 4);
				int legTwoStart = legOneEnd + BUFFER_REGION;
				legTwoStart = Math.min(legTwoStart, endIndex - 4);

				double sum = sliceLeg(index, bearings, times, MAX_ITERATIONS,
						legOneEnd, legTwoStart);

				thisSeries.add(new FixedMillisecond(times.get(index)), sum);
				straightBar.add(new FixedMillisecond(times.get(index)),
						overallScore);

				// is this better?
				if (sum < bestScore) {
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
		ChartPanel cp = new ChartPanel(new JFreeChart(
				"Results for " + SCENARIO, JFreeChart.DEFAULT_TITLE_FONT,
				combinedPlot, true)) {

			/**
					 * 
					 */
			private static final long serialVersionUID = 1L;

			@Override
			@Transient
			public Dimension getPreferredSize() {
				return new Dimension(1000, 800);
			}

		};
		frame.add(cp);
		frame.pack();

		long elapsed = System.currentTimeMillis() - startTime;
		System.out.println("Elapsed:" + elapsed / 1000 + " secs");

		// System.exit(0);
	}

	/**
	 * @return a frame to contain the results
	 */
	private static JFrame createFrame() {
		JFrame frame = new JFrame("Results");
		frame.pack();
		frame.setVisible(true);
		frame.addWindowListener(new WindowAdapter() {
			@Override
			public void windowClosing(WindowEvent e) {
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
	 * @param overallScore
	 *            the overall score for this leg
	 * @param BUFFER_REGION
	 * @param straightBar
	 * @param thisSeries
	 * @return
	 */
	private static double sliceLeg(int trialIndex, List<Double> bearings,
			List<Long> times, int MAX_ITERATIONS, final int legOneEnd,
			final int legTwoStart) {
		List<Long> theseTimes = times;
		List<Double> theseBearings = bearings;

		// first the times
		List<Long> beforeTimes = theseTimes.subList(0, legOneEnd);
		List<Long> afterTimes = theseTimes.subList(legTwoStart,
				theseTimes.size() - 1);

		List<Double> beforeBearings = theseBearings.subList(0, legOneEnd);
		List<Double> afterBearings = theseBearings.subList(legTwoStart,
				theseTimes.size() - 1);

		// fit a curve to the period after the turn
		PointValuePair beforeOptimiser = optimiseThis(beforeTimes,
				beforeBearings, beforeBearings.get(0));
		PointValuePair afterOptimiser = optimiseThis(afterTimes, afterBearings,
				afterBearings.get(0));

		// find the total error sum
		double sum = beforeOptimiser.getValue() / beforeTimes.size()  + afterOptimiser.getValue() / afterTimes.size();

		DecimalFormat intF = new DecimalFormat("00");
		
		if (trialIndex > 20 && trialIndex < 50) {
			System.out.println("index:" + intF.format(trialIndex)
//					+ " time:" + times.get(trialIndex) 
					+ " " + " Sum:" + numF.format(sum)
					+ " index:" + dateF.format(new Date(times.get(trialIndex)))
					+ " before:" + outDates(beforeTimes) + out(beforeOptimiser)
					+ " num:" + intF.format(beforeTimes.size())
			//		+ " after:" + outDates(afterTimes) 
			//		+ out(afterOptimiser)
			//		+ " num:" + intF.format(afterTimes.size())
);
		}

		return sum;
	}

	static PointValuePair optimiseThis(List<Long> times, List<Double> bearings,
			double initialBearing) {
		final MultivariateFunction function = new ArcTanSolver(times, bearings);
	
//		final SimplexOptimizer optimizerMult = new SimplexOptimizer(1e-4, 1e-7);
//		final SimplexOptimizer optimizerMult = new SimplexOptimizer(1e-3, 1e-6);
//		final SimplexOptimizer optimizerMult = new SimplexOptimizer(.0001,.0001);

//		final SimplexOptimizer optimizerMult = new SimplexOptimizer(1e-1, 1e-3);		
		final SimplexOptimizer optimizerMult = new SimplexOptimizer(1e-3, 1e-6);
		
		return optimizerMult.optimize(new MaxEval(MAX_ITERATIONS),
				new ObjectiveFunction(function), GoalType.MINIMIZE,
				new InitialGuess(new double[] { initialBearing, 0, 0 }),// afterBearings.get(0)
				new MultiDirectionalSimplex(3));
	}

	private static String outDates(List<Long> times) {
		String res = dateF.format(times.get(0)) + "-"
				+ dateF.format(times.get(times.size() - 1));
		return res;
	}

	//
	// System.out.println(name + " times:" + outDates(times) + " B:" +
	// numF.format(key[0]) + " P:" + numF.format(key[1]) + " Q:"
	// + numF.format(key[2]) + " Sum:" + wholeLegOptimiser.getValue());

	public static String out(PointValuePair res) {
		double[] key = res.getKey();
		String out = " B:" + numF.format(key[0]) + " P:"
				+ numF.format(key[1]) + " Q:"
				+ numF.format(key[2]) + " Sum:"
				+ numF.format(res.getValue());

		return out;
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
	private static List<LegOfData> calculateLegs(Track track) {

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

		for (int i = 0; i < times.length; i++) {
			long thisTime = times[i];

			double thisSpeed = speeds[i];
			double thisCourse = courses[i];

			if (i > 0) {
				// ok, check out the course change rate
				double timeStepSecs = (thisTime - lastTime) / 1000;
				double courseRate = Math.abs(thisCourse - lastCourse)
						/ timeStepSecs;
				double speedRate = Math.abs(thisSpeed - lastSpeed)
						/ timeStepSecs;

				// are they out of range
				if ((courseRate < COURSE_TOLERANCE)
						&& (speedRate < SPEED_TOLERANCE)) {
					// ok, we're on a new leg - drop the current one
					legs.get(legs.size() - 1).add(thisTime);
				} else {
					// we may be in a turn. create a new leg, if we haven't done
					// so already
					if (legs.get(legs.size() - 1).size() != 0) {
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

	/**
	 * function to generate sum of squares of errors for a single permutation of
	 * B,P,Q
	 * 
	 */
	private static class ArcTanSolver implements MultivariateFunction {
		final private List<Long> _times;
		final private List<Double> _bearings;

		public ArcTanSolver(List<Long> beforeTimes, List<Double> beforeBearings) {
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
				long elapsedMillis = _times.get(i) - _times.get(0);
				double elapsedSecs = elapsedMillis / 1000d;
				double thisForecast = calcForecast(B, P, Q, elapsedSecs);
				double thisMeasured = _bearings.get(i);
				double thisError = Math.pow(thisForecast - thisMeasured, 2);
				runningSum += thisError;

			}

			// System.out.println("B:" + (int)B + " P:" + P + " Q:" + Q +
			// " sum:" + runningSum);

			return runningSum / _times.size();
		}

		private double calcForecast(double B, double P, double Q,
				double elapsedSecs) {
			// return
			// Math.toDegrees(Math.atan2(Math.sin(Math.toRadians(B))+P*elapsedSecs,Math.cos(Math.toRadians(B))+Q*elapsedSecs));

			double dX = Math.cos(Math.toRadians(B)) + Q * elapsedSecs;
			double dY = Math.sin(Math.toRadians(B)) + P * elapsedSecs;

			// System.out.println("DX:" + dX + " dY:" + dY);

			return Math.toDegrees(Math.atan2(dY, dX));
		}
	}

}
