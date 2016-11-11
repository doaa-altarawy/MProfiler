package clustering.dataTypes;

import java.util.Vector;

public class Tree 
{	
    protected Double data;
    protected Tree left,right;    
    private Vector levelOrder;
    public int xl, xr, y;

    public Tree(){
        data = null;
        left = right = null;        
        levelOrder = new Vector(0,1);
    }
    public Tree(Double d){
        data = d;
        left = right = null;       
        levelOrder = new Vector(0,1);
    }
    public Tree(Double d, Tree left, Tree right){
        data = d;
        this.left = left;
        this.right = right;       
        levelOrder = new Vector(0,1);
    }
  
    public Vector getLevelOrder()
    {    	  	
    	getLevelOrder(this);
       	
    	
    	int n = levelOrder.size();
    	Tree temp;
    	for (int i=0; i<n/2+1; i++)
    	{
    		temp = (Tree)levelOrder.get(n-i-1);
    		levelOrder.set(n-i-1, levelOrder.get(i));
    		levelOrder.set(i, temp);
    	}
    	
    	return levelOrder;
    }
    
    private void getLevelOrder(Tree tree)
    {
    	levelOrder.add(tree);
    	if (tree.left!=null) levelOrder.add(tree.left);
    	if (tree.right!=null) levelOrder.add(tree.right);
    	
    	if (tree.left!=null) getLevelOrder(tree.left);
    	if (tree.right!=null) getLevelOrder(tree.right);
    	
    	    	
    }
    
    public static double getMaxNodeValue(Tree tree, double max)
    {
    	
    	if (tree.data!=null && tree.data>max) 
    		max = tree.data;
    	
    	if (tree.left==null && tree.right==null)
    		return max;
    	
    	double max1 = getMaxNodeValue(tree.left,max);
    	double max2 = getMaxNodeValue(tree.right,max);
    	if (max1>max2) return max1;
    	else return max2;
    	
    }
    
    public boolean isLeaf()
    {
    	if (left==null && right==null)
    		return true;
    	else return false;
    }
    
    public void setLeft(Tree l){
        left = l;
    }
    public void setRight(Tree r){
        right = r;
    }
    public void setData(Double d){
        data = d;
    }
    public Tree getLeft(){
        return left;
    }
    public Tree getRight(){
        return right;
    }
    public Double getData(){
        return data;
    }
		  

}
