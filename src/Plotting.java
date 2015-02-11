import java.awt.Color;
import java.text.SimpleDateFormat;

import javax.swing.JPanel;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.DateAxis;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.time.FixedMillisecond;
import org.jfree.data.time.TimeSeries;
import org.jfree.data.time.TimeSeriesCollection;

public class Plotting {

	public static void addOwnshipData(JPanel stack, Track ownshipTrack) {

		TimeSeriesCollection dataset1 = new TimeSeriesCollection();
		TimeSeriesCollection dataset2 = new TimeSeriesCollection();

		TimeSeries data1 = new TimeSeries("Course", FixedMillisecond.class);
		TimeSeries data2 = new TimeSeries("Speed", FixedMillisecond.class);

		double[] courses = ownshipTrack.getCourses();
		double[] speeds = ownshipTrack.getSpeeds();
		long[] times = ownshipTrack.getDates();

		// obtain the data for the points
		for (int i = 0; i < times.length; i++) {
			data1.add(new FixedMillisecond(times[i]), courses[i]);
			data2.add(new FixedMillisecond(times[i]), speeds[i]);
		}
		dataset1.addSeries(data1);
		dataset2.addSeries(data2);

		final JFreeChart chart = ChartFactory.createTimeSeriesChart				
				("Ownship", // String title,
				"Time", // String timeAxisLabel
				"Course", // String valueAxisLabel,
				dataset1, // XYDataset dataset,
				true, // include legend
				true, // tooltips
				false // urls
				);
		
		XYPlot xyPlot = (XYPlot) chart.getPlot();
		xyPlot.setDomainCrosshairVisible(true);
		xyPlot.setRangeCrosshairVisible(true);
	    final DateAxis axis = (DateAxis) xyPlot.getDomainAxis();
	    axis.setDateFormatOverride(new SimpleDateFormat("hh:mm:ss"));

	    final NumberAxis axis2 = new NumberAxis("Speed");
        xyPlot.setRangeAxis(1, axis2);
        xyPlot.setDataset(1, dataset2);
        xyPlot.mapDatasetToRangeAxis(1, 1);
	    
        XYLineAndShapeRenderer lineRenderer1 = new XYLineAndShapeRenderer(true, true);
        lineRenderer1.setSeriesPaint(1, Color.green);
        XYLineAndShapeRenderer lineRenderer2 = new XYLineAndShapeRenderer(true, true);
        lineRenderer2.setSeriesPaint(1, Color.blue);
        xyPlot.setRenderer(0, lineRenderer1);
        xyPlot.setRenderer(1, lineRenderer2);


        
//		xyPlot.setDataset(1, dataset2);
//
//		final NumberAxis xAxis2 = new NumberAxis("Speed");
//		xyPlot.setRangeAxis(1, xAxis2 );
////		xyPlot.mapDatasetToRangeAxis(0, 0);
//		xyPlot.mapDatasetToRangeAxis(1, 1);

		ChartPanel cp = new ChartPanel(chart);

		stack.add(cp);

	}

}
