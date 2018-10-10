# EMA
Electromagnetic Articulography Toolbox

HeadCorrection.R script will correct for head movement according to stable sensors on the sella-nasion plane. Sensors at the bridge of the nose (A), left mastoid (B), and right mastoid (C) create a plane at each frame of time. Vectors vBC, vMA (where M is the projection of A to line BC) and the cross product of vBC and vMA create new axes at each frame of time. Origin is set to the upper incisor (maxillary) sensor and primary articulators are redefined according to these axes.


TTPlotting.R uses ggplot to visualise the inferior-superior movement of the tongue tip for a given time frame, defined using the x-axis limits. A red horizontal line is super-imposed at 0, where the upper incisor sensor is located. This approximates raising to a critical "successful closure" point for a coronal stop.
