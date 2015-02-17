import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;


/** class to store a leg of ownship data
 * 
 * @author ian
 *
 */
public class LegOfData
{
	final List<Long> _times = new ArrayList<Long>();
	
	final private String _myName;
	
	public LegOfData(final String name)
	{
		_myName = name;
	}
	public Long getEnd() {
		return _times.get(_times.size()-1);
	}
	public Long getStart() {
		return _times.get(0);
	}
	public String getName() {
		return _myName;
	}
	public void add(long time)
	{
		_times.add(time);
	}
	public int size()
	{
		return _times.size();
	}
	@Override
	public String toString()
	{
		SimpleDateFormat sdf = new SimpleDateFormat("hh:mm:ss");
		return getName() + " " + sdf.format(new Date(_times.get(0))) + "-" + sdf.format(new Date(_times.get(_times.size()-1)));
	}
	
}