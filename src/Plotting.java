import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.text.SimpleDateFormat;
import java.util.Iterator;
import java.util.List;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.DateAxis;
import org.jfree.chart.axis.LogarithmicAxis;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.CombinedDomainXYPlot;
import org.jfree.chart.plot.IntervalMarker;
import org.jfree.chart.plot.Marker;
import org.jfree.chart.plot.ValueMarker;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.time.FixedMillisecond;
import org.jfree.data.time.TimeSeries;
import org.jfree.data.time.TimeSeriesCollection;
import org.jfree.ui.Layer;
import org.jfree.ui.RectangleAnchor;
import org.jfree.ui.TextAnchor;
import org.jfree.util.ShapeUtilities;

public class Plotting {

	public static CombinedDomainXYPlot createPlot()
	{
		final CombinedDomainXYPlot plot = new CombinedDomainXYPlot(new DateAxis("Domain"));
        plot.setGap(10.0);
		final DateAxis axis = (DateAxis) plot.getDomainAxis();
		axis.setDateFormatOverride(new SimpleDateFormat("hh:mm:ss"));
        return plot;
	}
	
	public static void addOwnshipData( CombinedDomainXYPlot parent, String title,
			Track ownshipTrack, List<LegOfData> ownshipLegs, 
			Color ownshipCol, 
			Track tgtTrack, List<LegOfData> tgtLegs, 
			Color tgtCol, 
			List<Long> turnEstimates, 
			final Long endTime) {

		TimeSeriesCollection courseColl = new TimeSeriesCollection();
		TimeSeriesCollection speedColl = new TimeSeriesCollection();

		TimeSeries oCourse = new TimeSeries("O/S Course", FixedMillisecond.class);
		TimeSeries oSpeed = new TimeSeries("O/S Speed", FixedMillisecond.class);
		TimeSeries tCourse = new TimeSeries("Tgt Course", FixedMillisecond.class);
		TimeSeries tSpeed = new TimeSeries("Tgt Speed", FixedMillisecond.class);

		double[] oCourses = ownshipTrack.getCourses();
		double[] oSpeeds = ownshipTrack.getSpeeds();
		long[] oTimes = ownshipTrack.getDates();
		double[] tCourses = tgtTrack.getCourses();
		double[] tSpeeds = tgtTrack.getSpeeds();

		// obtain the data for the points
		for (int i = 0; i < oTimes.length; i++) {
			long thisTime = oTimes[i];
			if(endTime == null || thisTime <= endTime)
			{
				oCourse.add(new FixedMillisecond(thisTime), oCourses[i]);
				oSpeed.add(new FixedMillisecond(thisTime), oSpeeds[i]);
				tCourse.add(new FixedMillisecond(thisTime), tCourses[i]);
				tSpeed.add(new FixedMillisecond(thisTime), tSpeeds[i]);
			}
		}
		courseColl.addSeries(oCourse);
		speedColl.addSeries(oSpeed);
		courseColl.addSeries(tCourse);
		speedColl.addSeries(tSpeed);

		final JFreeChart chart = ChartFactory.createTimeSeriesChart(title, // String
				"Time", // String timeAxisLabel
				title + "Course", // String valueAxisLabel,
				courseColl, // XYDataset dataset,
				true, // include legend
				true, // tooltips
				false); // urls 

		XYPlot xyPlot = (XYPlot) chart.getPlot();
		xyPlot.setDomainCrosshairVisible(true);
		xyPlot.setRangeCrosshairVisible(true);
		final DateAxis axis = (DateAxis) xyPlot.getDomainAxis();
		axis.setDateFormatOverride(new SimpleDateFormat("hh:mm:ss"));

		final NumberAxis axis2 = new NumberAxis(title + "Speed");
		xyPlot.setRangeAxis(1, axis2);
		xyPlot.setDataset(1, speedColl);
		xyPlot.mapDatasetToRangeAxis(1, 1);

		XYLineAndShapeRenderer lineRenderer1 = new XYLineAndShapeRenderer(true,
				true);
		lineRenderer1.setSeriesPaint(1, Color.blue);
		lineRenderer1.setSeriesShape(1, ShapeUtilities.createDiagonalCross(2, 2));
		lineRenderer1.setSeriesPaint(2, Color.red);
		lineRenderer1.setSeriesPaint(3, Color.blue);
		lineRenderer1.setSeriesPaint(4, Color.blue);
		lineRenderer1.setSeriesShape(4, ShapeUtilities.createDiagonalCross(2, 2));
		xyPlot.setRenderer(0, lineRenderer1);
		

		// let's try the shading
		if (ownshipLegs != null) {
			Iterator<LegOfData> iter = ownshipLegs.iterator();
			while (iter.hasNext()) {
				LegOfData leg = (LegOfData) iter.next();
				final Marker bst = new IntervalMarker(leg.getStart(),
						leg.getEnd(), ownshipCol, new BasicStroke(2.0f), null, null,
						1.0f);
				bst.setLabel(leg.getName());
				bst.setLabelAnchor(RectangleAnchor.BOTTOM_RIGHT);
				bst.setLabelFont(new Font("SansSerif", Font.ITALIC + Font.BOLD,
						10));
				bst.setLabelTextAnchor(TextAnchor.BASELINE_RIGHT);
				xyPlot.addDomainMarker(bst, Layer.BACKGROUND);
			}
		}
		
		if (tgtLegs != null) {
			Iterator<LegOfData> iter = tgtLegs.iterator();
			while (iter.hasNext()) {
				LegOfData leg = (LegOfData) iter.next();
				final Marker bst = new IntervalMarker(leg.getStart(),
						leg.getEnd(), tgtCol, new BasicStroke(2.0f), null, null,
						1.0f);
				bst.setLabel(leg.getName());
				bst.setLabelAnchor(RectangleAnchor.BOTTOM_RIGHT);
				bst.setLabelFont(new Font("SansSerif", Font.ITALIC + Font.BOLD,
						10));
				bst.setLabelTextAnchor(TextAnchor.BASELINE_RIGHT);
				xyPlot.addDomainMarker(bst, Layer.BACKGROUND);
			}
		}
		
		
		
		// let's try the shading
		if (turnEstimates != null) {
			plotMarkers(xyPlot, turnEstimates);
		}
		
		parent.add(xyPlot);
	}

	public static void addLegResults( CombinedDomainXYPlot parent, TimeSeriesCollection errorValues, List<Long> valueMarkers) {

		final JFreeChart chart = ChartFactory.createTimeSeriesChart(
				"Leg Results", // String
								// title,
				"Time", // String timeAxisLabel
				"Errpr", // String valueAxisLabel,
				errorValues, // XYDataset dataset,
				true, // include legend
				true, // tooltips
				false); // urls

		XYPlot xyPlot = (XYPlot) chart.getPlot();
		xyPlot.setDomainCrosshairVisible(true);
		xyPlot.setRangeCrosshairVisible(true);
		final DateAxis axis = (DateAxis) xyPlot.getDomainAxis();
		axis.setDateFormatOverride(new SimpleDateFormat("hh:mm:ss"));

		final NumberAxis rangeAxis = new LogarithmicAxis("Log(error)");
		xyPlot.setRangeAxis(rangeAxis);

		XYLineAndShapeRenderer lineRenderer1 = new XYLineAndShapeRenderer(true,
				true);
		xyPlot.setRenderer(0, lineRenderer1);
		
		// let's try the shading
		if (valueMarkers != null) {
			plotMarkers(xyPlot, valueMarkers);
		}
		
		xyPlot.getRenderer().setSeriesVisibleInLegend(false);
		
		
		parent.add(xyPlot);
	}

	/**  Plot a series of vertical markers
	 * @param xyPlot
	 * @param valueMarkers
	 */
	private static void plotMarkers(XYPlot xyPlot, List<Long> valueMarkers) {
		Iterator<Long> iter = valueMarkers.iterator();
		while (iter.hasNext()) {
			Long leg = (Long) iter.next();
			final Marker bst = new ValueMarker(leg,
					Color.gray, new BasicStroke(3.0f), null, null,
					1.0f);
			bst.setLabelAnchor(RectangleAnchor.BOTTOM_RIGHT);
			bst.setLabelFont(new Font("SansSerif", Font.ITALIC + Font.BOLD,
					10));
			bst.setLabelTextAnchor(TextAnchor.BASELINE_RIGHT);
			xyPlot.addDomainMarker(bst, Layer.BACKGROUND);
		}
	}
}
