package org.biojava.bio.tools;

/* not implemented */

public class MST
{
	private final int	maxN = 30;	// the max number of terminals
	private int		n = 10;		// the current number of terminals
	private final int	r = 4;		// the radius of a terminal
	private Point		p[];		// the terminals
	private Point		current;	// the current one and its old location
	private boolean		m[][];		// the minimum spanning tree edges
	
	private void mst() {
		  int dist[], neigh[], closest, minDist, d;

		  dist = new int[n];
		  neigh = new int[n];

		  // initialize data structures
		  for (int i = 0; i < n; i++) {
		    dist[i] = distance(p[0].x, p[0].y, p[i].x, p[i].y);
		    neigh[i] = 0;
		    for (int j = 0; j < n; j++) {
		      m[i][j] = false;
		    }
		  }

		  // find terminal closest to current partial tree
		  for (int i = 1; i < n; i++) {
		    closest = -1;
		    minDist = Integer.MAX_VALUE;
		    for (int j = 1; j < n; j++) {
		      if ((dist[j] != 0) && (dist[j] < minDist)) {
			closest = j;
			minDist = dist[j];
		      }
		    }
		    
		    // set an edge from it to its nearest neighbor
		    m[neigh[closest]][closest] = true;
		    m[closest][neigh[closest]] = true;

		    // update nearest distances to current partial tree
		    for (int j = 1; j < n; j++) {
		      d = distance(p[j].x, p[j].y, p[closest].x, p[closest].y);
		      if (d < dist[j]) {
			dist[j] = d;
			neigh[j] = closest;
		      }
		    }
		  }
		}
	
	private int distance(double x, double y, double x2, double y2)
	{
		// TODO Auto-generated method stub
		return 0;
	}

		// mst()
	
	private class Point{
		double x,y;
	}

}


