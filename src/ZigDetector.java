import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;

import junit.framework.TestCase;
import flanagan.math.Minimisation;
import flanagan.math.MinimisationFunction;

public class ZigDetector
{

	static class FlanaganArctan implements MinimisationFunction
	{
		private static double calcForecast(final double B, final double P,
				final double Q, final double elapsedSecs)
		{
			final double dX = Math.cos(Math.toRadians(B)) + Q * elapsedSecs;
			final double dY = Math.sin(Math.toRadians(B)) + P * elapsedSecs;
			return Math.toDegrees(Math.atan2(dY, dX));
		}

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
		@Override
		public double function(final double[] point)
		{
			final double B = point[0];
			final double P = point[1];
			final double Q = point[2];

			double runningSum = 0;
			final Long firstTime = _times[0];

			// ok, loop through the data
			for (int i = 0; i < _times.length; i++)
			{
				final long elapsedMillis = _times[i] - firstTime;
				final double elapsedSecs = elapsedMillis / 1000d;
				final double thisForecast = calcForecast(B, P, Q, elapsedSecs);
				final double thisMeasured = _bearings[i];
				double thisError = thisForecast - thisMeasured;
				if (thisError > 180)
				{
					thisError -= 360;
				}
				else if (thisError < -180)
				{
					thisError += 360;
				}
				final double sqError = Math.pow(thisError, 2);
				runningSum += sqError;
			}
			final double mean = runningSum / _times.length;

			final double rms = Math.sqrt(mean);
			return rms;
		}

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

	private final static int BUFFER_SIZE = 5;

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

	final SimpleDateFormat dateF = new SimpleDateFormat("HH:mm:ss");

	public static int getEnd(final int start, final int end, final int buffer,
			final int index)
	{
		int res = -1;
		final int MIN_SIZE = 3;
		final int semiBuffer = buffer / 2;

		final int earliestPossibleStartPoint = start + (MIN_SIZE - 1);

		if (index - semiBuffer > earliestPossibleStartPoint)
		{
			res = index - semiBuffer - 1;
		}
		return res;
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

	public static int getStart(final int start, final int end, final int buffer,
			final int index)
	{
		int res = -1;
		final int MIN_SIZE = 3;
		final int semiBuffer = buffer / 2;

		final int lastPossibleStartPoint = end - (MIN_SIZE - 1);

		if (index + semiBuffer < lastPossibleStartPoint)
		{
			res = index + semiBuffer + 1;
		}
		return res;
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

	/**
	 * slice this data into ownship legs, where the course and speed are
	 * relatively steady
	 * 
	 * @param course_degs
	 * @param speed
	 * @param bearings
	 * @param elapsedTimess
	 * @return
	 */
	public List<LegOfData> identifyOwnshipLegs(final Track track,
			final int avgPeriod)
	{

		final double COURSE_TOLERANCE = 0.1; // degs / sec (just a guess!!)
		final double SPEED_TOLERANCE = 0.01; // knots / sec (just a guess!!)

		double lastCourse = 0;
		double lastSpeed = 0;
		long lastTime = 0;

		final List<LegOfData> legs = new ArrayList<LegOfData>();
		legs.add(new LegOfData("Leg-1"));

		final long[] times = track.getDates();
		double[] speeds = track.getSpeeds();
		double[] courses = track.getCourses();

		// switch the courses to an n-term moving average
		courses = movingAverage(courses, avgPeriod);
		track.averageCourses = courses;

		speeds = movingAverage(speeds, avgPeriod);
		track.averageSpeeds = speeds;

		for (int i = 0; i < times.length; i++)
		{
			final long thisTime = times[i];

			final double thisSpeed = speeds[i];
			final double thisCourse = courses[i];

			if (i > 0)
			{
				// ok, check out the course change rate
				final double timeStepSecs = (thisTime - lastTime) / 1000;
				final double courseRate = Math.abs(thisCourse - lastCourse)
						/ timeStepSecs;
				final double speedRate = Math.abs(thisSpeed - lastSpeed) / timeStepSecs;

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
	private double[] movingAverage(final double[] measurements,
			final int period)
	{
		final double[] res = new double[measurements.length];
		final CenteredMovingAverage ma = new CenteredMovingAverage(period);
		for (int j = 0; j < measurements.length; j++)
		{
			// double d = measurements[j];
			// ma.newNum(d);
			// res[j] = ma.getAvg();
			res[j] = ma.average(j, measurements);
		}
		return res;
	}

	final private Minimisation optimiseThis(final List<Long> times,
			final List<Double> bearings, final double initialBearing,
			final double optimiserTolerance)
	{
		// Create instance of Minimisation
		final Minimisation min = new Minimisation();

		// Create instace of class holding function to be minimised
		final FlanaganArctan funct = new FlanaganArctan(times, bearings);

		// initial estimates
		final double firstBearing = bearings.get(0);
		final double[] start =
		{ firstBearing, 0.0D, 0.0D };

		// initial step sizes
		final double[] step =
		{ 0.2D, 0.3D, 0.3D };

		// convergence tolerance
		final double ftol = optimiserTolerance;

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
	 * @param optimiserTolerance
	 * @param fittedQ
	 * @param fittedP
	 * @param overallScore
	 *          the overall score for this leg
	 * @param BUFFER_REGION
	 * @param straightBar
	 * @param thisSeries
	 * @return
	 */
	private double sliceLeg(final int trialIndex,
			final List<Double> bearings, final List<Long> times,
			final int bufferSize, final int legOneEnd, final int legTwoStart,
			final double optimiserTolerance)
	{

		final List<Long> theseTimes = times;
		final List<Double> theseBearings = bearings;

		final Date thisD = new Date(times.get(trialIndex));

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
			final List<Long> beforeTimes = theseTimes.subList(0, legOneEnd);
			final List<Double> beforeBearings = theseBearings.subList(0, legOneEnd);
			beforeOptimiser = optimiseThis(beforeTimes, beforeBearings,
					beforeBearings.get(0), optimiserTolerance);
			beforeScore = beforeOptimiser.getMinimum();
			msg += " BEFORE:" + dateF.format(times.get(0)) + "-"
					+ dateF.format(times.get(legOneEnd)) + " " + beforeScore;
			// System.out.println(" before:" + _outDates(beforeTimes));
		}

		if (legTwoStart != -1)
		{
			final List<Long> afterTimes = theseTimes.subList(legTwoStart,
					theseTimes.size() - 1);
			final List<Double> afterBearings = theseBearings.subList(legTwoStart,
					theseTimes.size() - 1);
			afterOptimiser = optimiseThis(afterTimes, afterBearings,
					afterBearings.get(0), optimiserTolerance);
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
			final int beforeLen = theseTimes.subList(0, legOneEnd).size();
			final int afterLen = theseTimes.subList(legTwoStart,
					theseTimes.size() - 1).size();

			final int totalCuts = beforeLen + afterLen;

			final double beforeNormal = beforeScore * beforeLen / totalCuts;
			final double afterNormal = afterScore * afterLen / totalCuts;
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

	public void sliceThis(final String scenario, final long wholeStart,
			final long wholeEnd, final Sensor sensor, final ILegStorer legStorer,
			final double RMS_ZIG_RATIO, final double optimiseTolerance)
	{
		System.out.println("  Trying to slice : "
				+ dateF.format(new Date(wholeStart)) + " - "
				+ dateF.format(new Date(wholeEnd)));// + " " + curStart + " to " +
		// curEnd);

		// ok, find the best slice
		// prepare the data
		final List<Double> thisLegBearings = sensor.extractBearings(wholeStart,
				wholeEnd);
		final List<Long> thisLegTimes = sensor.extractTimes(wholeStart, wholeEnd);

		if (thisLegBearings.size() == 0)
		{
			return;
		}

		final Minimisation wholeLeg = optimiseThis(thisLegTimes, thisLegBearings,
				thisLegBearings.get(0), optimiseTolerance);
		final double wholeLegScore = wholeLeg.getMinimum();

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
			final int legOneEnd = getEnd(0, thisLegTimes.size(), BUFFER_SIZE, index);
			final int legTwoStart = getStart(0, thisLegTimes.size(), BUFFER_SIZE,
					index);

			// what's the total score for slicing at this index?
			final double sum = sliceLeg(index, thisLegBearings, thisLegTimes,
					BUFFER_SIZE, legOneEnd, legTwoStart, optimiseTolerance);

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
						bestScore / wholeLegScore * 100);
				legStorer.storeLeg(scenario, bestLegTwoStart, wholeEnd, sensor,
						bestScore / wholeLegScore * 100);
			}
			else
			{
				// right - we couldn't get a good slice. see what the whole score is
				System.out.println("Couldn't slice: whole leg score:" + wholeLegScore
						+ " best slice:" + bestScore);

				// just store the whole slice
				legStorer.storeLeg(scenario, wholeStart, wholeEnd, sensor,
						wholeLegScore / wholeLegScore * 100);
			}
		}
		else
		{
			System.out.println("slicing complete, can't slice");
		}
	}

}
