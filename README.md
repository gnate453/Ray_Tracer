Ray_Tracer
==========

As of 11/8/14 this program appears to meet the a3 specifications.

fixed: 
	-ambient lighting.
	-distacne to object origins and intersections measured from world coordinate of pixels instead of prp.
	-clipping planes (see note Issue 1)


Issue 1: 
	Some strange linear algebra required to fix near clipping plane. If an object is beyond the near clip plane (Object(Z) < Pixel(Z)) then d^2, c^2, and v are calculated as lecture slides suggest ( (Origin - Pixel) and (r^2 - (c^2-v^2)) ), however if the object is in front of the near clipping plane (Object(Z) > Pixel (Z)) then d^2 c^2 and v are calculated in reverse ( (Pixel(Z) - Origin) and ((c^2 + v^2) - r^2) ). I'm not sure why this leads me to correct near clipping behavior.  I arrived at this solution by experimentation, visualually inspection each modifcation and adjusting to the corret behavior.

	

