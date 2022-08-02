package localsearch.domainspecific.vehiclerouting.apps.MTDLCVR.service;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;


public class SubSetGenerator {
	private Set<Integer> S;
	public ArrayList<Integer> map;
	private int n;
	private ArrayList<Integer> X;
	public HashSet<HashSet<Integer>> subSet;
	
	public SubSetGenerator(Set<Integer> S){
		this.S = S;
		map = new ArrayList<Integer>();
		X = new ArrayList<Integer>();
		subSet = new HashSet<HashSet<Integer>>();
		
		this.n = S.size();
		
		for(int e : S){
			map.add(e);
			X.add(0);
		}
		
		
	}
	
	public void solution(){
		HashSet<Integer> s = new HashSet<Integer>();
		for(int i = 0; i < X.size(); i++){
			if(X.get(i) == 1)
				s.add(map.get(i));
		}
		subSet.add(s);
	}
	
	public void TRY(int k){
		for(int v = 0; v < 2; v++){
			X.set(k, v);
			if(k == n - 1)
				solution();
			else
				TRY(k+1);
		}
	}
	
	public void generate(){
		TRY(0);
	}
}
