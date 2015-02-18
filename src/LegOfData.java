import java.text.SimpleDateFormat;
import java.util.Date;


/** class to store a leg of ownship data
 * 
 * @author ian
 *
 */
public class LegOfData
{
	long tStart = Long.MAX_VALUE;
	long tEnd = Long.MIN_VALUE;
	
	final private String _myName;
	
	public LegOfData(final String name)
	{
		_myName = name;
	}
	public LegOfData(String name, long tStart2, long tEnd2)
	{
		this(name);
		tStart = tStart2;
		tEnd = tEnd2;
	}
	public boolean initialised()
	{
		return tStart != Long.MAX_VALUE;
	}
	public Long getEnd() {
		return tEnd;
	}
	public Long getStart() {
		return tStart;
	}
	public String getName() {
		return _myName;
	}
	public void add(long time)
	{
		tStart = Math.min(tStart, time);
		tEnd = Math.max(tEnd, time);
	}
	@Override
	public String toString()
	{
		SimpleDateFormat sdf = new SimpleDateFormat("hh:mm:ss");
		return getName() + " " + sdf.format(new Date(tStart)) + "-" + sdf.format(new Date(tEnd));
	}
	
}