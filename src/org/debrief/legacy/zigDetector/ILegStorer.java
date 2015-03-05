package org.debrief.legacy.zigDetector;

public interface ILegStorer
{

	/**
	 * register this leg of data
	 * 
	 * @param scenario
	 * @param tStart
	 * @param tEnd
	 */
	void storeLeg(String scenarioName, long tStart, long tEnd, Sensor sensor,
			double rms);

}