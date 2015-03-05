# ZigDetector
Algorithm to detect target zigs

The following screenshot shows a summary of a single run of the algorithm.  The individual plots are described beneath the plot.
![Overview results](http://i.imgur.com/NL84H2e.png)

## Plot One
This is a plot of the course & speed of the ownship vessel. The green sectors mark the periods where ownship is on steady course & speed. They are labelled *Leg-1*, *Leg-2* and *Leg-3*.

## Plot Two
This is a plot of the error sums produced when experimenting with slicing the legs at different points. An line is also shown showing the error sum for the whole of that leg.

The vertical markers show the slice time that produces the minimum error score.

## Plot Three
This is a plot of the target course & speed. It is clear that the turn at around 12:37 has been identified correctly, but not the turn at around 01:22

The above strategy culminated with the multi-leg consumer algorithm stored in this commit: https://github.com/IanMayo/ZigDetector/commit/192bd79f4f2eda60b57ed4583e58ad8aaec8482b

In that strategy we walked along the target data, storing a leg each time we got a block of data that met an error threshold.

But, this strategy was only of use in nice, neat simulate data. It totally fell over once real, dirty data was obtained - since it produced lots of very small threshold-compliant legs.  For that, we had to return to the original ArcTan strategy - where we optimised slicing target legs by comparing the normalised sum of the before/after legs against the RMS score for the whole leg.

# Real data
Hey, how does that crazy algorithm perform with real-world, dirty, smelly data?

This next image shows the results when using very poor quality ownship course/speed data, with low fidelity quantised bearing data.

![Poor data](http://i.imgur.com/AUVtu0B.png)

There are some tune-able variables available that have been used to correctly slice the ownship legs.  But slicing the target legs pretty-much worked out of the box.
