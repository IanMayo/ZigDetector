import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.text.SimpleDateFormat;
import java.util.Iterator;
import java.util.List;

import javax.swing.JPanel;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.DateAxis;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.IntervalMarker;
import org.jfree.chart.plot.Marker;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.time.FixedMillisecond;
import org.jfree.data.time.TimeSeries;
import org.jfree.data.time.TimeSeriesCollection;
import org.jfree.ui.Layer;
import org.jfree.ui.RectangleAnchor;
import org.jfree.ui.TextAnchor;

public class Plotting {

	public static void addOwnshipData(JPanel stack, String title, Track ownshipTrack,
			List<LegOfData> ownshipLegs) {

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

		final JFreeChart chart = ChartFactory.createTimeSeriesChart(title, // String
																				// title,
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

		XYLineAndShapeRenderer lineRenderer1 = new XYLineAndShapeRenderer(true,
				true);
		lineRenderer1.setSeriesPaint(1, Color.green);
		XYLineAndShapeRenderer lineRenderer2 = new XYLineAndShapeRenderer(true,
				true);
		lineRenderer2.setSeriesPaint(1, Color.blue);
		xyPlot.setRenderer(0, lineRenderer1);
		xyPlot.setRenderer(1, lineRenderer2);

		// let's try the shading
		if (ownshipLegs != null) {
			Iterator<LegOfData> iter = ownshipLegs.iterator();
			while (iter.hasNext()) {
				LegOfData leg = (LegOfData) iter.next();
				final Color c = new Color(55, 255, 24, 63);
				final Marker bst = new IntervalMarker(leg.getStart(),
						leg.getEnd(), c, new BasicStroke(2.0f), null, null,
						1.0f);
				bst.setLabel(leg.getName());
				bst.setLabelAnchor(RectangleAnchor.BOTTOM_RIGHT);
				bst.setLabelFont(new Font("SansSerif", Font.ITALIC + Font.BOLD,
						10));
				bst.setLabelTextAnchor(TextAnchor.BASELINE_RIGHT);
				xyPlot.addDomainMarker(bst, Layer.BACKGROUND);
			}
		}

		ChartPanel cp = new ChartPanel(chart);
		stack.add(cp);
	}
}
