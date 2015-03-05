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

import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.CombinedDomainXYPlot;
import org.jfree.chart.plot.XYPlot;
import org.jfree.data.time.FixedMillisecond;
import org.jfree.data.time.TimeSeries;
import org.jfree.data.time.TimeSeriesCollection;

import flanagan.interpolation.CubicSpline;
import flanagan.math.Minimisation;

public class ZigDetectorTest
{

	private static final String OPTIMISER_THRESHOLD_STR = "Optim Thresh";

	private static final String RMS_ZIG_RATIO_STR = "RMS Zig Ratio";

	private static final String NOISE_SD_STR = "Noise SD";

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

//		panel.add(createItem(NOISE_SD_STR, new double[]
//		{ 0d, 0.1d, 0.2d, 0.25d, 0.3d, 0.5, 2d }, newListener, decFormat,
//				BRG_ERROR_SD));
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

		final ZigDetectorTest detector = new ZigDetectorTest();

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

						ZigDetector.sliceThis(data._name, legStart, legEnd, data.sensor,
								data.legStorer, RMS_ZIG_RATIO, OPTIMISER_TOLERANCE);

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
			data.ownshipLegs = ZigDetector.identifyOwnshipLegs(data.ownshipTrack, 5);
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

	/**
	 * local instance of leg storer, also collects some other performance data
	 * 
	 * @author ian
	 * 
	 */
	public static class LegStorer implements ILegStorer
	{
		private ArrayList<LegOfData> _legList;
		private TimeSeries _rmsScores;

		@Override
		public void storeLeg(String scenarioName, long tStart, long tEnd,
				Sensor sensor, double rms)
		{
			System.out.println("Storing " + scenarioName + " : "
					+ dateF.format(new Date(tStart)) + " - "
					+ dateF.format(new Date(tEnd)));

			_legList.add(new LegOfData("Leg-" + (_legList.size() + 1), tStart, tEnd));

			// store some RMS error scores
			if (sensor != null)
			{
				List<Long> times = sensor.extractTimes(tStart, tEnd);
				for (Iterator<Long> iterator = times.iterator(); iterator.hasNext();)
				{
					Long long1 = (Long) iterator.next();
					_rmsScores.add(new FixedMillisecond(long1), rms);
				}
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

		@Override
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

	// private static void sliceThisConsumer(String scenario, long curStart,
	// long curEnd, Sensor sensor, LegStorer legStorer)
	// {
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

	// }

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
