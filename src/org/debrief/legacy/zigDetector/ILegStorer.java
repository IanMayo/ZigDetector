package org.debrief.legacy.zigDetector;

/** interface for listener class that is told when a new leg is detected
 * 
 * @author ian
 *
 */
public interface ILegStorer
{

	/**
	 * register this leg of data
	 * 
	 * @param scenario
	 * @param tStart
	 * @param tEnd
	 * @param sensor - optional sensor object that produced the data
	 * @param rms - the %age error from the RMS for the whole leg
	 */
	void storeLeg(String scenarioName, long tStart, long tEnd, Sensor sensor,
			double rms);

}