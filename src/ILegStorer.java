import java.util.List;


public interface ILegStorer
{

	/** register this leg of data
	 * 
	 * @param scenario
	 * @param tStart
	 * @param tEnd
	 */
	void storeLeg(String scenarioName, long tStart, long tEnd, Sensor sensor,
			double rms);

	/** return the legs that have been found
	 * 
	 * @return
	 */
	List<LegOfData> getLegs();

	
}