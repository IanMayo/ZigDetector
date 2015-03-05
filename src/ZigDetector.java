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
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.beans.Transient;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;

import javax.swing.ButtonGroup;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JRadioButton;

import junit.framework.TestCase;

import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.CombinedDomainXYPlot;
import org.jfree.chart.plot.XYPlot;
import org.jfree.data.time.FixedMillisecond;
import org.jfree.data.time.TimeSeries;
import org.jfree.data.time.TimeSeriesCollection;

import flanagan.interpolation.CubicSpline;
import flanagan.math.Minimisation;
import flanagan.math.MinimisationFunction;

public class ZigDetector
{

	private static final String OPTIMISER_THRESHOLD_STR = "Optim Thresh";

	private static final String RMS_ZIG_RATIO_STR = "RMS Zig Ratio";

	private static final String NOISE_SD_STR = "Noise SD";

	private static final int BUFFER_SIZE = 5;

	// how much RMS error we require on the Atan Curve before we
	// bother trying to slice the target leg
	private static double RMS_ZIG_RATIO = 0.6;
	// private static double RMS_ZIG_THRESHOLD = 0.005;

	// when to let the optimiser relax
	private static double OPTIMISER_TOLERANCE = 1e-6;

	// the error we add on
	private static double BRG_ERROR_SD = 0.0;

	final static Long timeEnd = null; // osL1end;

	final static SimpleDateFormat dateF = new SimpleDateFormat("HH:mm:ss");
	final static DecimalFormat numF = new DecimalFormat(
			" 0000.0000000;-0000.0000000");

	protected JPanel createControls(NewValueListener newListener)
	{
		JPanel panel = new JPanel();
		panel.setLayout(new GridLayout(0, 1));

		NumberFormat decFormat = new DecimalFormat("0.000");
		NumberFormat expFormat = new DecimalFormat("0.000E00");

		panel.add(createItem(NOISE_SD_STR, new double[]
		{ 0d, 0.1d, 0.2d, 0.25d, 0.3d, 0.5, 2d }, newListener, decFormat,
				BRG_ERROR_SD));
		panel.add(createItem(OPTIMISER_THRESHOLD_STR, new double[]
		{ 1e-4, 1e-5, 1e-6, 1e-7, 1e-8 }, newListener, expFormat,
				OPTIMISER_TOLERANCE));
		panel.add(createItem(RMS_ZIG_RATIO_STR, new double[]
		{ 1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3 }, newListener, decFormat,
				RMS_ZIG_RATIO));

		return panel;
	}

	public static void main(String[] args) throws Exception
	{

		final ZigDetector detector = new ZigDetector();

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

		final HashMap<String, ScenDataset> datasets = new HashMap<String, ScenDataset>();

		ArrayList<String> scenarios = new ArrayList<String>();
		// scenarios.add("Scen1");
		// scenarios.add("Scen2a");
		// scenarios.add("Scen2b");
		// scenarios.add("Scen3");
		// scenarios.add("Scen4");
		scenarios.add("Scen5");

		// handler for slider changes
		NewValueListener newL = new NewValueListener()
		{

			@Override
			public void newValue(String name, double val)
			{
				switch (name)
				{
				case NOISE_SD_STR:
					BRG_ERROR_SD = val;
					break;
				case RMS_ZIG_RATIO_STR:
					RMS_ZIG_RATIO = val;
					break;
				case OPTIMISER_THRESHOLD_STR:
					OPTIMISER_TOLERANCE = val;
					break;
				default:
					// don't worry - we make it empty to force a refresh
				}

				// create the ownship & target course data
				Iterator<ScenDataset> iterator = datasets.values().iterator();
				while (iterator.hasNext())
				{
					// create somewhere to store it.
					final ScenDataset data = iterator.next();

					// clear the identified legs - to show progress
					Plotting.clearLegMarkers(data.targetPlot);

					// update the title
					NumberFormat numF = new DecimalFormat("0.000");
					NumberFormat expF = new DecimalFormat("0.###E0");
					String title = data._name + " Noise:" + numF.format(BRG_ERROR_SD)
							+ " Conv:" + expF.format(OPTIMISER_TOLERANCE) + " Zig:"
							+ numF.format(RMS_ZIG_RATIO);
					data.chartPanel.getChart().setTitle(title);

					data.turnMarkers = new ArrayList<Long>();

					// get ready to store the new legs
					data.legStorer = new LegStorer();

					// apply the error to the sensor data
					data.sensor.applyError(BRG_ERROR_SD);

					// get ready to store the results runs
					TimeSeriesCollection legResults = new TimeSeriesCollection();

					TimeSeries rmsScores = new TimeSeries("RMS Errors");
					data.legStorer.setRMSScores(rmsScores);
					data.legStorer.setLegList(new ArrayList<LegOfData>());

					// ok, work through the legs. In the absence of a Discrete
					// Optimisation
					// algorithm we're taking a brue force approach.
					// Hopefully Craig can find an optimised alternative to this.
					for (Iterator<LegOfData> iterator2 = data.ownshipLegs.iterator(); iterator2
							.hasNext();)
					{
						LegOfData thisLeg = (LegOfData) iterator2.next();

						// ok, slice the data for this leg
						long legStart = thisLeg.getStart();
						long legEnd = thisLeg.getEnd();

						// trim the start/end to the sensor data
						legStart = Math.max(legStart, data.sensor.getTimes()[0]);
						legEnd = Math.min(legEnd,
								data.sensor.getTimes()[data.sensor.getTimes().length - 1]);

						sliceThis(data._name, legStart, legEnd, data.sensor, data.legStorer);

						// create a placeholder for the overall score for this leg
						TimeSeries atanBar = new TimeSeries("ATan " + thisLeg.getName());
						legResults.addSeries(atanBar);

						// create a placeholder for the individual time slice experiments
						TimeSeries thisSeries = new TimeSeries(thisLeg.getName()
								+ " Slices");
						legResults.addSeries(thisSeries);
					}

					// plot the bearings
					Plotting.showBearings(data.bearingPlot, data.sensor.getTimes(),
							data.sensor.getBearings(), data.legStorer._rmsScores);

					// ok, output the results
					Plotting.plotLegPeriods(data.targetPlot, data.tgtTransColor,
							data.legStorer._legList);
				}
			}
		};

		// - ok insert the grid controls
		inGrid.add(detector.createControls(newL));

		// create the placeholders
		for (Iterator<String> iterator = scenarios.iterator(); iterator.hasNext();)
		{
			String name = (String) iterator.next();

			// create somewhere to store it.
			ScenDataset data = new ScenDataset(name);

			// store it
			datasets.put(name, data);

			// ok - create the placeholder
			CombinedDomainXYPlot combinedPlot = Plotting.createPlot();
			data._plot = combinedPlot;

			data.chartPanel = new ChartPanel(new JFreeChart("Results for " + name
					+ " Tol:" + OPTIMISER_TOLERANCE, JFreeChart.DEFAULT_TITLE_FONT,
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
			inGrid.add(data.chartPanel, BorderLayout.CENTER);
		}

		// create the ownship & target course data
		Iterator<ScenDataset> iterator = datasets.values().iterator();
		while (iterator.hasNext())
		{
			// create somewhere to store it.
			ScenDataset data = iterator.next();

			// load the data
			data.ownshipTrack = Track.read("data/" + data._name + "_Ownship.csv");
			data.targetTrack = Track.read("data/" + data._name + "_Target.csv");
			data.sensor = new Sensor("data/" + data._name + "_Sensor.csv");

			// find the ownship legs
			data.ownshipLegs = identifyOwnshipLegs(data.ownshipTrack, 5);
			// data.ownshipLegs = data.ownshipLegs.subList(2, 3);

			// ok, now for the ownship data
			data.oShipColor = new Color(0f, 0f, 1.0f);
			data.oShipTransColor = new Color(0f, 0f, 1.0f, 0.2f);
			data.ownshipPlot = Plotting.plotSingleVesselData(data._plot, "O/S",
					data.ownshipTrack, data.oShipColor, null, timeEnd);

			// ok, now for the ownship legs
			Plotting.plotLegPeriods(data.ownshipPlot, data.oShipTransColor,
					data.ownshipLegs);

			// try to plot the moving average
			// switch the courses to an n-term moving average
			Plotting.addAverageCourse(data.ownshipPlot,
					data.ownshipTrack.averageCourses, data.ownshipTrack.averageSpeeds,
					data.ownshipTrack.getDates());

			// and the target plot
			data.tgtColor = new Color(1.0f, 0f, 0f);
			data.tgtTransColor = new Color(1.0f, 0f, 0f, 0.2f);
			data.targetPlot = Plotting.plotSingleVesselData(data._plot, "Tgt",
					data.targetTrack, data.tgtColor, null, timeEnd);

			// insert a bearing plot
			data.bearingPlot = Plotting.createBearingPlot(data._plot);
		}

		if (inGrid.getComponentCount() == 1)
			grid.setColumns(1);

		frame.pack();

		// ok, we should probably initialise it
		newL.newValue("", 0);

		long elapsed = System.currentTimeMillis() - startTime;
		System.out.println("Elapsed:" + elapsed / 1000 + " secs");

	}

	public static class ScenDataset
	{
		public XYPlot bearingPlot;
		public ChartPanel chartPanel;
		public Color tgtTransColor;
		public Color tgtColor;
		public Color oShipColor;
		public Color oShipTransColor;
		protected ArrayList<Long> turnMarkers;
		protected LegStorer legStorer;
		public XYPlot targetPlot;
		public XYPlot ownshipPlot;
		public List<LegOfData> ownshipLegs;
		private String _name;
		CombinedDomainXYPlot _plot;
		public Track ownshipTrack;
		public Track targetTrack;
		public Sensor sensor;

		public ScenDataset(String name)
		{
			_name = name;
		}

	}

	public static class LegStorer
	{
		private ArrayList<LegOfData> _legList;
		private TimeSeries _rmsScores;

		public void storeLeg(String scenario, long tStart, long tEnd,
				Sensor sensor, double rms)
		{
			System.out.println("Storing " + scenario + " : "
					+ dateF.format(new Date(tStart)) + " - "
					+ dateF.format(new Date(tEnd)));

			_legList.add(new LegOfData("Leg-" + (_legList.size() + 1), tStart, tEnd));

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

	protected static interface ValConverter
	{
		public double convert(int input);

		int unConvert(double val);
	}

	protected static interface NewValueListener
	{
		public void newValue(String attribute, double val);
	}

	protected JPanel createItem(final String label, final double[] values,
			final NewValueListener listener, final NumberFormat numberFormat,
			double startValue)
	{

		JPanel panel = new JPanel();
		panel.setLayout(new BorderLayout());
		panel.add(new JLabel(label), BorderLayout.WEST);

		JPanel buttons = new JPanel();
		buttons.setLayout(new GridLayout(1, 0));
		ButtonGroup bg = new ButtonGroup();
		for (int i = 0; i < values.length; i++)
		{
			final JRadioButton newB = new JRadioButton();
			final double thisD = values[i];
			newB.setText(numberFormat.format(values[i]));
			if (thisD == startValue)
				newB.setSelected(true);
			newB.addActionListener(new ActionListener()
			{

				@Override
				public void actionPerformed(ActionEvent e)
				{
					System.out.println(label + " changed to:" + e.toString());

					if (newB.isSelected())
					{
						listener.newValue(label, thisD);
					}
				}
			});
			buttons.add(newB);
			bg.add(newB);
		}

		panel.add(buttons, BorderLayout.EAST);

		return panel;
	}

	private static void sliceThis(String scenario, final long wholeStart,
			final long wholeEnd, Sensor sensor, LegStorer legStorer)
	{
		System.out.println("  Trying to slice : "
				+ dateF.format(new Date(wholeStart)) + " - "
				+ dateF.format(new Date(wholeEnd)));// + " " + curStart + " to " +
		// curEnd);

		// ok, find the best slice
		// prepare the data
		List<Double> thisLegBearings = sensor.extractBearings(wholeStart, wholeEnd);
		List<Long> thisLegTimes = sensor.extractTimes(wholeStart, wholeEnd);

		if (thisLegBearings.size() == 0)
			return;

		Minimisation wholeLeg = optimiseThis(thisLegTimes, thisLegBearings,
				thisLegBearings.get(0));
		double wholeLegScore = wholeLeg.getMinimum();

		// ok, now have to slice it
		double bestScore = Double.MAX_VALUE;
		int bestSlice = -1;
		long sliceTime = -1;
		long bestLegOneEnd = -1;
		long bestLegTwoStart = -1;

		// find the optimal first slice
		for (int index = 0; index < thisLegTimes.size(); index++)
		// for (int index = 12; index < 17; index++)
		{
			int legOneEnd = getEnd(0, thisLegTimes.size(), BUFFER_SIZE, index);
			int legTwoStart = getStart(0, thisLegTimes.size(), BUFFER_SIZE, index);

			// what's the total score for slicing at this index?
			double sum = sliceLeg(index, thisLegBearings, thisLegTimes, BUFFER_SIZE,
					legOneEnd, legTwoStart);

			// is this better?
			if ((sum > 0) && (sum < bestScore))
			{
				// yes - store it.
				bestScore = sum;
				bestSlice = index;
				sliceTime = thisLegTimes.get(index);
				bestLegOneEnd = thisLegTimes.get(legOneEnd);
				bestLegTwoStart = thisLegTimes.get(legTwoStart);
			}
		}

		// right, how did we get on?
		if (sliceTime != -1)
		{
			System.out.println("  Best slice at:" + dateF.format(new Date(sliceTime))
					+ " index:" + bestSlice + " score:" + bestScore + " whole leg:"
					+ wholeLegScore);

			// is this slice acceptable?
			if (bestScore < wholeLegScore * RMS_ZIG_RATIO)
			{
				legStorer.storeLeg(scenario, wholeStart, bestLegOneEnd, sensor,
						bestScore);
				legStorer.storeLeg(scenario, bestLegTwoStart, wholeEnd, sensor,
						bestScore);
			}
			else
			{
				// right - we couldn't get a good slice. see what the whole score is
				System.out.println("Couldn't slice: whole leg score:" + wholeLegScore
						+ " best slice:" + bestScore);

				// just store the whole slice
				legStorer.storeLeg(scenario, wholeStart, wholeEnd, sensor,
						wholeLegScore);
			}
		}
		else
		{
			System.out.println("slicing complete, can't slice");
		}
	}

	private static void sliceThisConsumer(String scenario, long curStart,
			long curEnd, Sensor sensor, LegStorer legStorer)
	{
		// System.out.println("  Trying to slice : "
		// + dateF.format(new Date(curStart)) + " - "
		// + dateF.format(new Date(curEnd)));// + " " + curStart + " to " +
		// // curEnd);
		//
		// // ok, find the best slice
		// // prepare the data
		// List<Double> thisLegBearings = sensor.extractBearings(curStart, curEnd);
		// List<Long> thisLegTimes = sensor.extractTimes(curStart, curEnd);
		//
		// if (thisLegBearings.size() == 0)
		// return;
		//
		// Minimisation wholeLeg = optimiseThis(thisLegTimes, thisLegBearings,
		// thisLegBearings.get(0));
		// // is this slice acceptable?
		// if (wholeLeg.getMinimum() < RMS_ZIG_THRESHOLD)
		// {
		// legStorer.storeLeg(scenario, curStart, curEnd, sensor,
		// wholeLeg.getMinimum());
		// }
		// else
		// {
		// // ok, we'll have to slice it
		// double bestScore = Double.MAX_VALUE;
		// int bestSlice = -1;
		// long sliceTime = -1;
		//
		// // find the optimal first slice
		// for (int index = 0; index < thisLegTimes.size(); index++)
		// // for (int index = 12; index < 17; index++)
		// {
		// // what's the total score for slicing at this index?
		// double sum = sliceLeg(index, thisLegBearings, thisLegTimes, BUFFER_SIZE);
		//
		// // System.out.println("score at: " + dateF.format(new
		// Date(thisLegTimes.get(index))) + " = " + sum);
		//
		// // System.out.println("score for:" + dateF.format(new
		// // Date(thisLegTimes.get(index))) + " is " + sum);
		//
		// // is this better?
		// if ((sum > 0) && (sum < bestScore))
		// {
		// // yes - store it.
		// bestScore = sum;
		// bestSlice = index;
		// sliceTime = thisLegTimes.get(index);
		// }
		// }
		//
		// // right, how did we get on?
		// if (sliceTime != -1)
		// {
		// System.out.println("  Best slice at:"
		// + dateF.format(new Date(sliceTime)) + " index:" + bestSlice
		// + " score:" + bestScore);
		//
		// // is this slice acceptable?
		// if (bestScore < RMS_ZIG_THRESHOLD)
		// {
		// legStorer.storeLeg(scenario, curStart, sliceTime, sensor, bestScore);
		//
		// // move the leg start along a little, to allow for a turn
		// long newLegStart = sliceTime + 60000;
		//
		// // do we have enough to look at?
		// long remainingTime = curEnd - newLegStart;
		//
		// System.out.println("  about to slice, have " + (remainingTime / 1000) +
		// " secs left");
		//
		// // try to slice it
		// sliceThis(scenario, newLegStart, curEnd, sensor, legStorer);
		// }
		// else
		// {
		// List<Double> trimLegBearings = sensor.extractBearings(curStart,
		// sliceTime);
		// List<Long> trimLegTimes = sensor.extractTimes(curStart, sliceTime);
		//
		// // ok, see if we can reduce the buffer size
		// boolean found = false;
		// int bufferLen = 1;
		// int maxBuffer = Math.min(trimLegTimes.size() - 3, 7);
		//
		// while (!found && bufferLen < maxBuffer)
		// {
		// List<Long> beforeTimes = trimLegTimes.subList(0,
		// trimLegTimes.size() - bufferLen);
		// List<Double> beforeBearings = trimLegBearings.subList(0,
		// trimLegTimes.size() - bufferLen);
		//
		// // System.out.println(" trimming:" + _outDates(beforeTimes));
		//
		// Minimisation beforeOptimiser = optimiseThis(beforeTimes,
		// beforeBearings, beforeBearings.get(0));
		// double sum = beforeOptimiser.getMinimum();
		//
		// if (sum < RMS_ZIG_THRESHOLD)
		// {
		// // ok, we can move on.
		// found = true;
		// legStorer.storeLeg(scenario, curStart,
		// beforeTimes.get(beforeTimes.size() - 1), sensor, sum);
		// curStart = sliceTime + 1;
		// }
		//
		// bufferLen++;
		// }
		//
		// // ok, did we find anything?
		// if (found)
		// {
		// // ok, we'll automatically move along
		// sliceThis(scenario, curStart, curEnd, sensor, legStorer);
		// }
		// else
		// {
		// // right - it's damn impossible! force it to move along
		//
		// }
		// }
		// }
		// else
		// {
		// System.out.println("slicing complete, can't slice");
		// }
		// }

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
		double ftol = OPTIMISER_TOLERANCE;

		// Nelder and Mead minimisation procedure
		min.nelderMead(funct, start, step, ftol);

		return min;
	}

	/**
	 * @param trialIndex
	 * @param bearings
	 * @param times
	 * @param legOneEnd
	 * @param legTwoStart
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
			List<Long> times, int bufferSize, int legOneEnd, int legTwoStart)
	{

		List<Long> theseTimes = times;
		List<Double> theseBearings = bearings;

		Date thisD = new Date(times.get(trialIndex));

		// if((legOneEnd == -1) || (legTwoStart == -1))
		// return Double.MAX_VALUE;

		double beforeScore = Double.MAX_VALUE;
		double afterScore = Double.MAX_VALUE;

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
					+ dateF.format(times.get(legOneEnd)) + " " + beforeScore;
			// System.out.println(" before:" + _outDates(beforeTimes));
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
					+ dateF.format(times.get(times.size() - 1)) + " " + afterScore;
			// System.out.println(" after:" + _outDates(afterTimes));
		}

		// find the total error sum
		double sum = beforeScore + afterScore;

		// do we have both legs?
		if ((legOneEnd != -1) && (legTwoStart != -1))
		{
			int beforeLen = theseTimes.subList(0, legOneEnd).size();
			int afterLen = theseTimes.subList(legTwoStart, theseTimes.size() - 1)
					.size();

			int totalCuts = beforeLen + afterLen;

			double beforeNormal = beforeScore * beforeLen / totalCuts;
			double afterNormal = afterScore * afterLen / totalCuts;
			sum = beforeNormal + afterNormal;
		}

		// System.out.println(msg+ " SUM:" + sum);

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
	private static List<LegOfData> identifyOwnshipLegs(Track track, int avgPeriod)
	{

		final double COURSE_TOLERANCE = 0.1; // degs / sec (just a guess!!)
		final double SPEED_TOLERANCE = 0.01; // knots / sec (just a guess!!)

		double lastCourse = 0;
		double lastSpeed = 0;
		long lastTime = 0;

		List<LegOfData> legs = new ArrayList<LegOfData>();
		legs.add(new LegOfData("Leg-1"));

		long[] times = track.getDates();
		double[] speeds = track.getSpeeds();
		double[] courses = track.getCourses();

		// switch the courses to an n-term moving average
		courses = movingAverage(courses, avgPeriod);
		track.averageCourses = courses;

		speeds = movingAverage(speeds, avgPeriod);
		track.averageSpeeds = speeds;

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

	/**
	 * create a moving average over the set of dat measurements
	 * 
	 * @param measurements
	 * @param period
	 * @return
	 */
	private static double[] movingAverage(double[] measurements, int period)
	{
		double[] res = new double[measurements.length];
		CenteredMovingAverage ma = new CenteredMovingAverage(period);
		for (int j = 0; j < measurements.length; j++)
		{
			// double d = measurements[j];
			// ma.newNum(d);
			// res[j] = ma.getAvg();
			res[j] = ma.average(j, measurements);
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

	final private static class FlanaganArctan implements MinimisationFunction
	{
		final private Long[] _times;
		final private Double[] _bearings;

		public FlanaganArctan(final List<Long> beforeTimes,
				final List<Double> beforeBearings)
		{
			_times = beforeTimes.toArray(new Long[]
			{});
			_bearings = beforeBearings.toArray(new Double[]
			{});
		}

		// evaluation function
		public double function(double[] point)
		{
			final double B = point[0];
			final double P = point[1];
			final double Q = point[2];

			double runningSum = 0;
			final Long firstTime = _times[0];

			// ok, loop through the data
			for (int i = 0; i < _times.length; i++)
			{
				long elapsedMillis = _times[i] - firstTime;
				double elapsedSecs = elapsedMillis / 1000d;
				double thisForecast = calcForecast(B, P, Q, elapsedSecs);
				double thisMeasured = _bearings[i];
				double thisError = thisForecast - thisMeasured;
				if (thisError > 180)
				{
					thisError -= 360;
				}
				else if (thisError < -180)
				{
					thisError += 360;
				}
				double sqError = Math.pow(thisError, 2);
				runningSum += sqError;
			}
			double mean = runningSum / _times.length;

			double rms = Math.sqrt(mean);
			return rms;
		}

	}

	@SuppressWarnings("unused")
	private static double _splineErrorfor(List<Long> times,
			List<Double> bearings, Minimisation wholeLegOptimiser)
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
		TimeSeries p = new TimeSeries("P-calc");
		TimeSeries q = new TimeSeries("Q-calc");
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

	public static class MovingAverage
	{
		private final Queue<Double> window = new LinkedList<Double>();
		private final int period;
		private double sum;

		public MovingAverage(int period)
		{
			assert period > 0 : "Period must be a positive integer";
			this.period = period;
		}

		public void newNum(double num)
		{
			sum += num;
			window.add(num);
			if (window.size() > period)
			{
				sum -= window.remove();
			}
		}

		public double getAvg()
		{
			if (window.isEmpty())
				return 0; // technically the average is undefined
			return sum / window.size();
		}

	}
}
