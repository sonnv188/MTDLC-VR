package localsearch.domainspecific.vehiclerouting.apps.MTDLCVR.service;

import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Scanner;
import java.util.Set;

import gurobi.*;


public class MILPsolver{
	Set<Integer> truckPoints;
	
	Set<Integer> customerPoints;
	Set<Integer> centralDepots;
	Set<Integer> visitedPoints;
	
	int nbProducts;
	double[] weightProducts;
	
	int nbCustomers;
	double[][] demands;
	int[] earliestArrivalTime;
	int[] latestArrivalTime;
	int[] waitingDuration;
	double[] durationTimeUnit;
	double[] serviceDuration;
	
	int nbParkings;
	int[] nbVehiclesAtParking;
	ArrayList<Integer> parkings;
	HashMap<Integer, Integer> vh2parking;
	
	int nbCentralDepots;
	//ArrayList<Integer> depots;
	
	int nbVehicles;
	int[] vhStartWorkingTime;
	int[] vhEndWorkingTime;
	double[] vhLowerCapacity;
	double[] vhUpperCapacity;
	int[] vhNbTours;
	
	int[][] vhRestrictProducts;
	int[][] vhRestrictCustomers;
	int[][] vhRemainCustomers;
	
	double[][] travelTime;	
	
	ArrayList<GRBVar> X;//xe k di tu i  den j trong chuyen q
	HashMap<String, GRBVar> arc2X;
	
//	ArrayList<GRBVar> Y;//xe k phuc vu kh i trong chuyen q
//	HashMap<String, GRBVar> arc2Y;
	
	ArrayList<GRBVar> L;//so luong hang p load len xe k trong chuyen q tai depot
	HashMap<String, GRBVar> arc2L;
	
	
	ArrayList<GRBVar> Z;
	HashMap<String, GRBVar> arc2Z;
	
	ArrayList<GRBVar> W;
	HashMap<String, GRBVar> point2W;
	
	ArrayList<GRBVar> arrivalTime;
	ArrayList<GRBVar> departureTime;
	ArrayList<GRBVar> startServiceTime;
	HashMap<String, GRBVar> point2ArrivalTime;
	HashMap<String, GRBVar> point2DepartureTime;
	HashMap<String, GRBVar> point2StartServiceTime;
	
	ArrayList<GRBVar> servingTimeAtPoint;
	HashMap<String, GRBVar> point2servingTimeAtPoint;
	
	ArrayList<GRBVar> arrivalTimeDepot;
	ArrayList<GRBVar> returnTimeDepot;
	ArrayList<GRBVar> departureTimeDepot;
	ArrayList<GRBVar> startServiceTimeDepot;
	HashMap<String, GRBVar> point2ArrivalTimeDepot;
	HashMap<String, GRBVar> point2returnTimeDepot;
	HashMap<String, GRBVar> point2DepartureTimeDepot;
	HashMap<String, GRBVar> point2StartServiceTimeDepot;
	
	HashSet<ArrayList<Integer>> arcVehicles;
	HashMap<Integer, HashSet<Integer>> inArcVehicles;
	HashMap<Integer, HashSet<Integer>> outArcVehicles;
	
	GRBVar y0;
	GRBVar y1;
	
	int a = 100;
	int M = 57600;
	int M1 = 57600;
	int M2 = 57600;
	int M3 = 57600;
	int M4 = 57600;
	int M5 = 57600;
	int M6 = 57600;
	int M7 = 57600;
	
	public MILPsolver(){
	}
	
	public void readData(String fn){
		truckPoints = new HashSet<Integer>();
		visitedPoints = new HashSet<Integer>();
		customerPoints = new HashSet<Integer>();
		centralDepots = new HashSet<Integer>();
		parkings = new ArrayList<Integer>();
		
		int id = 0;
		
		try{
			Scanner sc = new Scanner(new File(fn));
			while(sc.hasNextLine()){
				String[] str;
				System.out.println(sc.nextLine());
				nbCustomers = Integer.parseInt(sc.nextLine());
				
				sc.nextLine();
				nbParkings = Integer.parseInt(sc.nextLine());
				
				sc.nextLine();
				nbCentralDepots = Integer.parseInt(sc.nextLine());
				
				sc.nextLine();
				nbVehicles = Integer.parseInt(sc.nextLine());
				
				sc.nextLine();
				nbProducts = Integer.parseInt(sc.nextLine());
				
				sc.nextLine();
				nbVehiclesAtParking	= new int[nbParkings];
				earliestArrivalTime = new int[nbCustomers+nbParkings+nbCentralDepots];
				latestArrivalTime 	= new int[nbCustomers+nbParkings+nbCentralDepots];
				waitingDuration		= new int[nbCustomers+nbParkings+nbCentralDepots];
				durationTimeUnit 	= new double[nbCustomers+nbParkings+nbCentralDepots];
				
				for(int i = 0; i < nbParkings; i++){
					str = sc.nextLine().split(" ");
					nbVehiclesAtParking[id] = Integer.parseInt(str[0]);
					earliestArrivalTime[id] = Integer.parseInt(str[1]);
					latestArrivalTime[id] = Integer.parseInt(str[2]);
					waitingDuration[id] = 0;
					durationTimeUnit[id] = 0;
					truckPoints.add(id);
					parkings.add(id);
					id++;
				}
				
				sc.nextLine();
				
				
				//depots = new ArrayList<Integer>();
				for(int i = 0; i < nbCentralDepots; i++){
					str = sc.nextLine().split(" ");
					truckPoints.add(id);
					//depots.add(id);
					centralDepots.add(id);
					visitedPoints.add(id);
					earliestArrivalTime[id] = Integer.parseInt(str[0]);
					latestArrivalTime[id] = Integer.parseInt(str[1]);
					waitingDuration[id] = Integer.parseInt(str[2]);
					durationTimeUnit[id] = Double.parseDouble(str[3]);
					id++;
				}
				
				sc.nextLine();
				vhStartWorkingTime = new int[nbVehicles];
				vhEndWorkingTime = new int[nbVehicles];
				vhLowerCapacity = new double[nbVehicles];
				vhUpperCapacity = new double[nbVehicles];
				vhNbTours = new int[nbVehicles];
				vh2parking = new HashMap<Integer, Integer>();
				
				int k = 0;
				for(int i = 0; i < nbParkings; i++){
					for(int j = 0; j < nbVehiclesAtParking[i]; j++){
						str = sc.nextLine().split(" ");
						vhStartWorkingTime[k] = Integer.parseInt(str[0]);
						vhEndWorkingTime[k] = Integer.parseInt(str[1]);
						vhLowerCapacity[k] = Double.parseDouble(str[2]);
						vhUpperCapacity[k] = Double.parseDouble(str[3]);
						vhNbTours[k] = Integer.parseInt(str[4]);
						vh2parking.put(k, i);
						k++;
					}
				}
			
				sc.nextLine();
				weightProducts = new double[nbProducts];
				for(int i = 0; i < nbProducts; i++)
					weightProducts[i] = Double.parseDouble(sc.nextLine());
				
				sc.nextLine();
				demands = new double[nbCustomers+nbParkings+nbCentralDepots][nbProducts];
				for(int i = 0; i < nbParkings; i++)
					for(int j = 0; j < nbProducts; j++)
						demands[i][j] = 0;
				for(int i = 0; i < nbCustomers; i++){
					str = sc.nextLine().split(" ");
					for(int j = 0; j < nbProducts; j++)
						demands[i+nbParkings+nbCentralDepots][j] = Double.parseDouble(str[j]);
				}
				
				sc.nextLine();
				
				for(int i = 0; i < nbCustomers; i++){
					str = sc.nextLine().split(" ");
					truckPoints.add(id);
					customerPoints.add(id);
					visitedPoints.add(id);
					earliestArrivalTime[id] = Integer.parseInt(str[0]);
					latestArrivalTime[id] = Integer.parseInt(str[1]);
					waitingDuration[id] = Integer.parseInt(str[2]);
					durationTimeUnit[id] = Double.parseDouble(str[3]);
					id++;
				}
				
				sc.nextLine();
				vhRestrictProducts = new int[nbVehicles][nbProducts];
				for(int i = 0; i < nbVehicles; i++){
					str = sc.nextLine().split(" ");
					for(int j = 0; j < nbProducts; j++)
						vhRestrictProducts[i][j] = Integer.parseInt(str[j]);
				}
				
				sc.nextLine();
				vhRestrictCustomers = new int[nbVehicles][nbCustomers+nbParkings+nbCentralDepots];
				for(int i = 0; i < nbVehicles; i++)
					for(int j = 0; j < nbParkings + nbCentralDepots; j++)
						vhRestrictCustomers[i][j] = 1;
				for(int i = 0; i < nbVehicles; i++){
					str = sc.nextLine().split(" ");
					for(int j = 0; j < nbCustomers; j++)
						vhRestrictCustomers[i][j+nbParkings+nbCentralDepots] = Integer.parseInt(str[j]);
				}
				
				sc.nextLine();
				vhRemainCustomers = new int[nbVehicles][nbCustomers+nbParkings+nbCentralDepots];
				for(int i = 0; i < nbVehicles; i++)
					for(int j = 0; j < nbParkings; j++)
						vhRemainCustomers[i][j] = 0;
				for(int i = 0; i < nbVehicles; i++){
					str = sc.nextLine().split(" ");
					for(int j = 0; j < nbCustomers; j++)
						vhRemainCustomers[i][j+nbParkings+nbCentralDepots] = Integer.parseInt(str[j]);
				}
				
				sc.nextLine();
				travelTime = new double[truckPoints.size()][truckPoints.size()];
				for(int i = 0; i < truckPoints.size(); i++)
					for(int j = 0; j < truckPoints.size(); j++)
						travelTime[i][j] = 0;
				int nbLines = Integer.parseInt(sc.nextLine());
				for(int i = 0; i < nbLines; i++){
					str = sc.nextLine().split(" ");
					int f = Integer.parseInt(str[0]);
					int t = Integer.parseInt(str[1]);
					double d = Double.parseDouble(str[2]);
					travelTime[f][t] = d;
				}
				
			}
			sc.close();
		}catch(Exception e){
			System.out.println(e);
		}
		
		createArcTrucks();
		calculateServeDurationAtCustomer();
	}
	
	public void createArcTrucks(){
		arcVehicles = new HashSet<ArrayList<Integer>>();
		inArcVehicles = new HashMap<Integer, HashSet<Integer>>();
		outArcVehicles = new HashMap<Integer, HashSet<Integer>>();
		
		for(int u : truckPoints){
			for(int v : truckPoints){
				if(u == v)
					continue;
				if(parkings.contains(u) && parkings.contains(v))
					continue;
				if(parkings.contains(u) && customerPoints.contains(v))
					continue;
				if(centralDepots.contains(u) && parkings.contains(v))
					continue;
				if(centralDepots.contains(u) && centralDepots.contains(v))
					continue;
//				if(customerPoints.contains(u) && centralDepots.contains(v))
//					continue;
				ArrayList<Integer> arc = new ArrayList<Integer>();
				arc.add(u);
				arc.add(v);
				arcVehicles.add(arc);
			}
		}
		
		for(int v : truckPoints){
			inArcVehicles.put(v, new HashSet<Integer>());
			outArcVehicles.put(v, new HashSet<Integer>());
		}

		for(ArrayList<Integer> arc : arcVehicles){
			HashSet<Integer> outKey = outArcVehicles.get(arc.get(0));
			outKey.add(arc.get(1));
			outArcVehicles.put(arc.get(0), outKey);
			HashSet<Integer> inVal = inArcVehicles.get(arc.get(1));
			inVal.add(arc.get(0));
			inArcVehicles.put(arc.get(1), inVal);
		}
	}
	
	public void defineFlowVariableOfVehicles(GRBModel model){
		X = new ArrayList<GRBVar>();
		arc2X = new HashMap<String, GRBVar>();
		
		for(int k = 0; k < nbVehicles; k++){
			for(int q = 0; q < vhNbTours[k]; q++){
				for(ArrayList<Integer> arc : arcVehicles){
					int u = arc.get(0);
					int v = arc.get(1);
					try {
						GRBVar var = model.addVar(0.0, 1.0, 0.0, GRB.BINARY, "X(" + k + "," + u + "," + v + "," + q + ")");
						X.add(var);
						String s = k + "-" + u + "-" + v + "-" + q;
						arc2X.put(s, var);
					} catch (GRBException e) {
						System.out.println("Error code: " + e.getErrorCode() + ". " +
					            e.getMessage());
					}
				}
			}
		}
	}
	
	public void defineQuantityOfLoadAtDepot(GRBModel model){
		L = new ArrayList<GRBVar>();
		arc2L = new HashMap<String, GRBVar>();
		
		for(int k = 0; k < nbVehicles; k++){
			for(int q = 0; q < vhNbTours[k]; q++){
				for(int p = 0; p < nbProducts; p++){
					int maxItem = (int)(vhUpperCapacity[k]/weightProducts[p]);
					try {
						GRBVar var = model.addVar(0.0, maxItem, 0.0, GRB.INTEGER, "L(" + k + "," + p + "," + q + ")");
						L.add(var);
						String s = k + "-" + p + "-" + q;
						arc2L.put(s, var);
					} catch (GRBException e) {
						System.out.println("Error code: " + e.getErrorCode() + ". " +
					            e.getMessage());
					}
				}
			}
		}
	}
	
	public void defineOperationOfTrip(GRBModel model){
		Z = new ArrayList<GRBVar>();
		arc2Z = new HashMap<String, GRBVar>();
		
		for(int k = 0; k < nbVehicles; k++){
			for(int q = 0; q < vhNbTours[k]; q++){
				try {
					GRBVar var = model.addVar(0.0, 1.0, 0.0, GRB.BINARY, "Z(" + k + "," + q + ")");
					Z.add(var);
					String s = k + "-" + q;
					arc2Z.put(s, var);
				} catch (GRBException e) {
					System.out.println("Error code: " + e.getErrorCode() + ". " +
				            e.getMessage());
				}
			}
		}
	}
	
	public void defineLastPointOfTripVariables(GRBModel model){
		W = new ArrayList<GRBVar>();
		point2W = new HashMap<String, GRBVar>();
		
		for(int k = 0; k < nbVehicles; k++){
			for(int q = 0; q < vhNbTours[k]; q++){
				for(int i : customerPoints){
					try {
						GRBVar var = model.addVar(0.0, 1.0, 0.0, GRB.BINARY, "W(" + k + "," + i + "," + q + ")");
						W.add(var);
						String s = k + "-" + i + "-" + q;
						point2W.put(s, var);
					} catch (GRBException e) {
						System.out.println("Error code: " + e.getErrorCode() + ". " +
					            e.getMessage());
					}
				}
			}
		}
	}
	
	public void defineServingTimeAtPointVariables(GRBModel model){
		servingTimeAtPoint = new ArrayList<GRBVar>();
		point2servingTimeAtPoint = new HashMap<String, GRBVar>();
		
		for(int k = 0; k < nbVehicles; k++){
			for(int i : truckPoints){
				for(int q = 0; q < vhNbTours[k]; q++){
					try{
						GRBVar var = model.addVar(0.0, vhEndWorkingTime[k], 0.0, GRB.INTEGER, "servingTimeAtPoint(" + k + "," + i + "," + q + ")");
						servingTimeAtPoint.add(var);
						String key = k + "-" + i + "-" + q;
						point2servingTimeAtPoint.put(key, var);
					} catch (GRBException e) {
						System.out.println("Error code: " + e.getErrorCode() + ". " +
					            e.getMessage());
					}
				}
			}
		}
	}
	
	public void defineBinaryVariables(GRBModel model){
		try{
			y0 = model.addVar(0.0, 1, 0.0, GRB.BINARY, "y0");			
			y1 = model.addVar(0.0, 1, 0.0, GRB.BINARY, "y1");
			
		} catch (GRBException e) {
			System.out.println("Error code: " + e.getErrorCode() + ". " +
		            e.getMessage());
		}
	}
	
	public void defineVariables(GRBModel model){
		defineFlowVariableOfVehicles(model);
//		defineServedAtCustomerOfVehicles(model);
		defineQuantityOfLoadAtDepot(model);
		defineOperationOfTrip(model);
		defineLastPointOfTripVariables(model);
//		defineArrivalDepartureTimeVariables(model);
//		defineStartServiceTimeVariables(model);
//		defineArrivalDepartureTimeDepotVariables(model);
//		defineStartServiceTimeDepotVariables(model);
//		defineReturnTimeDepotVariables(model);
		defineBinaryVariables(model);
		defineServingTimeAtPointVariables(model);
	}
	
	public void flowBalanceConstraint(GRBModel model){
		try{
			for(int k = 0; k < nbVehicles; k++){
				for(int u : truckPoints){
					GRBLinExpr expr = new GRBLinExpr();
					for(int q = 0; q < vhNbTours[k]; q++){
						for(int v : inArcVehicles.get(u)){
							String s = k + "-" + v + "-" + u + "-" + q;
							GRBVar x = arc2X.get(s);
							expr.addTerm(1, x);
						}
						for(int v : outArcVehicles.get(u)){
							String s = k + "-" + u + "-" + v + "-" + q;
							GRBVar x = arc2X.get(s);
							expr.addTerm(-1, x);
						}
					}
					model.addConstr(expr, GRB.EQUAL, 0.0, "BalanceFlow(" + k + "," + u + ")");
				}
			}
		} catch (GRBException e) {
			System.out.println("Error code: " + e.getErrorCode() + ". " +
		            e.getMessage());
		}
	}
	
	public void flowBalanceConstraintMiller(GRBModel model){
		try{
			for(int u : customerPoints){
				GRBLinExpr expr = new GRBLinExpr();
				GRBLinExpr expr2 = new GRBLinExpr();
				for(int k = 0; k < nbVehicles; k++){
					for(int q = 0; q < vhNbTours[k]; q++){
						for(int v : inArcVehicles.get(u)){
							String s = k + "-" + v + "-" + u + "-" + q;
							GRBVar x = arc2X.get(s);
							expr.addTerm(1, x);
						}
						for(int v : outArcVehicles.get(u)){
							String s = k + "-" + u + "-" + v + "-" + q;
							GRBVar x = arc2X.get(s);
							expr2.addTerm(1, x);
						}
					}
				}
				model.addConstr(expr, GRB.EQUAL, 1.0, "BalanceFlow1(" + u + ")");
				model.addConstr(expr2, GRB.EQUAL, 1.0, "BalanceFlow2(" + u + ")");
			}
		} catch (GRBException e) {
			System.out.println("Error code: " + e.getErrorCode() + ". " +
		            e.getMessage());
		}
	}
	
	public void flowBalanceOnEachTripConstraint(GRBModel model){
		try{
			for(int k = 0; k < nbVehicles; k++){
				for(int q = 0; q < vhNbTours[k]; q++){
					for(int u : customerPoints){
						GRBLinExpr expr = new GRBLinExpr();
						for(int v : inArcVehicles.get(u)){
							String s = k + "-" + v + "-" + u + "-" + q;
							GRBVar x = arc2X.get(s);
							expr.addTerm(1, x);
						}
						for(int v : outArcVehicles.get(u)){
							String s = k + "-" + u + "-" + v + "-" + q;
							GRBVar x = arc2X.get(s);
							expr.addTerm(-1, x);
						}
						model.addConstr(expr, GRB.EQUAL, 0.0, "BalanceFlow(" + k + "," + u + ")");
					}
				}
			}
		} catch (GRBException e) {
			System.out.println("Error code: " + e.getErrorCode() + ". " +
		            e.getMessage());
		}
	}
	
	public void SubtourEliminationConstraint(GRBModel model){
		SubSetGenerator ssgen = new SubSetGenerator(visitedPoints);
		ssgen.generate();
		int t = 0;
		try{
			for(HashSet<Integer> s : ssgen.subSet){
				HashSet<ArrayList<Integer>> A = getArcInduced(s, arcVehicles);
				int m = s.size() - 1;
				if(m <= 0)
					continue;
				for(int k = 0; k < nbVehicles; k++){
					GRBLinExpr expr = new GRBLinExpr();	
					for(int q = 0; q < vhNbTours[k]; q++){	
						for(ArrayList<Integer> arc : arcVehicles){
							int au = arc.get(0);
							int av = arc.get(1);
							String str = k + "-" + au + "-" + av + "-" + q;
							if(arc2X.get(str) == null)
								continue;
							
							GRBVar x = arc2X.get(str);
							if(checkExist(A, arc))
								expr.addTerm(1, x);
						}
					}
					model.addConstr(expr, GRB.LESS_EQUAL, m, "SubTour(" + k + "," + t + ")");
					t++;
				}
			}
		} catch (GRBException e) {
			System.out.println("Error code: " + e.getErrorCode() + ". " +
		            e.getMessage());
		}
	}
	
	public void SubtourEliminationInTripConstraint(GRBModel model){
		SubSetGenerator ssgen = new SubSetGenerator(visitedPoints);
		ssgen.generate();
		int t = 0;
		try{
			for(HashSet<Integer> s : ssgen.subSet){
				HashSet<ArrayList<Integer>> A = getArcInduced(s, arcVehicles);
				int m = s.size() - 1;
				if(m <= 0)
					continue;
				for(int k = 0; k < nbVehicles; k++){
					for(int q = 0; q < vhNbTours[k]; q++){
						GRBLinExpr expr = new GRBLinExpr();	
						for(ArrayList<Integer> arc : arcVehicles){
							int au = arc.get(0);
							int av = arc.get(1);
							String str = k + "-" + au + "-" + av + "-" + q;
							if(arc2X.get(str) == null)
								continue;
							
							GRBVar x = arc2X.get(str);
							if(checkExist(A, arc))
								expr.addTerm(1, x);
						}
						model.addConstr(expr, GRB.LESS_EQUAL, m, "SubTourInTrip(" + k + "," + t + ")");
						t++;
					}
				}
			}
		} catch (GRBException e) {
			System.out.println("Error code: " + e.getErrorCode() + ". " +
		            e.getMessage());
		}
	}
	
	public void defineEachCustomerVisitedAtMostOnceConstraint(GRBModel model){
		try{
			for(int i : customerPoints){
				GRBLinExpr expr = new GRBLinExpr();
				for(int k = 0; k < nbVehicles; k++){
					for(int q = 0; q < vhNbTours[k]; q++){
						for(ArrayList<Integer> arc : arcVehicles){
							int au = arc.get(0);
							int av = arc.get(1);
							String s = k + "-" + au + "-" + av + "-" + q;
							GRBVar x = arc2X.get(s);
							if(av == i){
								expr.addTerm(1, x);
							}
						}
					}
				}
				model.addConstr(expr, GRB.EQUAL, 1.0, "CustomerVisitedOnce(" + i + ")");
			}
		} catch (GRBException e) {
			System.out.println("Error code: " + e.getErrorCode() + ". " +
		            e.getMessage());
		}
	}
	
	public void defineOneDepotVisitedOneTimeForEachTrip(GRBModel model){
		try{
			for(int k = 0; k < nbVehicles; k++){
				for(int q = 0; q < vhNbTours[k]; q++){
					GRBLinExpr expr1 = new GRBLinExpr();
					if(q == 0){
						for(int dp : centralDepots){
							for(int j : parkings){
								String kjdpq = k + "-" + j + "-" + dp + "-" + q;
								GRBVar x1 = arc2X.get(kjdpq);
								expr1.addTerm(1, x1);
							}
						}
						
						String kq = k + "-" + q;
						GRBVar z = arc2Z.get(kq);
						expr1.addTerm(-1, z);
						model.addConstr(expr1, GRB.EQUAL, 0.0, "oneDepotVisitedOnetimeForEachTrip(" + k + "," + q + ")");
					}
					else{
						for(int dp : centralDepots){
							for(int j : customerPoints){
								String kjdpq1 = k + "-" + j + "-" + dp + "-" + (q-1);
								GRBVar x1 = arc2X.get(kjdpq1);
								expr1.addTerm(1, x1);
							}
						}
						
						String kq = k + "-" + q;
						GRBVar z = arc2Z.get(kq);
						expr1.addTerm(-1, z);
						model.addConstr(expr1, GRB.EQUAL, 0.0, "oneDepotVisitedOnetimeForEachTrip(" + k + "," + q + ")");
					}
				}
			}
		} catch (GRBException e) {
			System.out.println("Error code: " + e.getErrorCode() + ". " +
		            e.getMessage());
		}
	}
	
	public void defineEachVehicleNotStartAtParkingForTheNextTrip(GRBModel model){
		try{
			for(int k = 0; k < nbVehicles; k++){
//				//start at parking 1 lan toi 1 central depot
//				GRBLinExpr expr1 = new GRBLinExpr();
//				
//				String k0 = k + "-0";
//				GRBVar z = arc2Z.get(k0);
//				expr1.addTerm(-1, z);
//				int pk = vh2parking.get(k);
//				for(int dp : centralDepots){
//					String kpkdp0 = k + "-" + pk + "-" + dp + "-0";
//					GRBVar x = arc2X.get(kpkdp0);
//					expr1.addTerm(1, x);
//				}
//				model.addConstr(expr1, GRB.EQUAL, 0.0, "vhStartAtParkingOneTime(" + k + ")");
				
				//khong start tai cac trip tiep theo
				GRBLinExpr expr2 = new GRBLinExpr();
				int pk = vh2parking.get(k);
				for(int q = 1; q < vhNbTours[k]; q++){
					for(int dp : centralDepots){
						String kpkdpq = k + "-" + pk + "-" + dp + "-" + q;
						GRBVar x = arc2X.get(kpkdpq);
						expr2.addTerm(1, x);
					}
				}
				model.addConstr(expr2, GRB.EQUAL, 0.0, "vhNotStartAtParkingForTheNextTrip(" + k + "," + pk + ")");
			}
		} catch (GRBException e) {
			System.out.println("Error code: " + e.getErrorCode() + ". " +
		            e.getMessage());
		}
	}
	
	public void defineVehicleNotStartEndAtOtherParkingConstraint(GRBModel model){		
		try{
			//not start
			for(int k = 0; k < nbVehicles; k++){
				for(int q = 0; q < vhNbTours[k]; q++){
					for(int pk : parkings){
						if(pk == vh2parking.get(k))
							continue;
						GRBLinExpr expr1 = new GRBLinExpr();
						for(int j : centralDepots){
							String k0jq = k + "-" + pk + "-" + j + "-" + q;
							GRBVar x1 = arc2X.get(k0jq);
							expr1.addTerm(1, x1);
						}
						model.addConstr(expr1, GRB.EQUAL, 0.0, "vhNotStartAtOtherParking1(" + k + "," + pk + "," + q + ")");
					}
				}
			}
			
			//not end
			for(int k = 0; k < nbVehicles; k++){
				for(int q = 0; q < vhNbTours[k]; q++){
					for(int pk : parkings){
						if(pk == vh2parking.get(k))
							continue;
						GRBLinExpr expr2 = new GRBLinExpr();
						for(int j : customerPoints){							
							String kj0q = k + "-" + j + "-" + pk + "-" + q;
							GRBVar x2 = arc2X.get(kj0q);
							expr2.addTerm(1, x2);
						}
						model.addConstr(expr2, GRB.EQUAL, 0.0, "vhNotStartAtOtherParking2(" + k + "," + pk + "," + q + ")");
					}
				}
			}
		} catch (GRBException e) {
			System.out.println("Error code: " + e.getErrorCode() + ". " +
		            e.getMessage());
		}
	}
	
	public void defineEachVehicleStartsAndComesBackParkingAtMostOnce(GRBModel model){
		try{
			for(int k = 0; k < nbVehicles; k++){
				GRBLinExpr expr1 = new GRBLinExpr();
				GRBLinExpr expr2 = new GRBLinExpr();
				int pk = vh2parking.get(k);
				for(int j : centralDepots){
					String kpkjq = k + "-" + pk + "-" + j + "-0";
					GRBVar x1 = arc2X.get(kpkjq);
					expr1.addTerm(1, x1);
				}
				String k0 = k + "-0";
				GRBVar z = arc2Z.get(k0);
				model.addConstr(expr1, GRB.EQUAL, z, "vhStartParkingAtMostone1(" + k + "," + pk + ")");
				//model.addConstr(expr2, GRB.EQUAL, z, "vhNotStartAtOtherParking2(" + k + "," + pk + ")");
			}
			for(int k = 0; k < nbVehicles; k++){
				GRBLinExpr expr2 = new GRBLinExpr();
				int pk = vh2parking.get(k);
				for(int q = 0; q < vhNbTours[k]; q++) {
					for(int j : customerPoints){
						String kjpkq = k + "-" + j + "-" + pk + "-" + q;
						GRBVar x1 = arc2X.get(kjpkq);
						expr2.addTerm(1, x1);
					}
				}
				String k0 = k + "-0";
				GRBVar z = arc2Z.get(k0);
				model.addConstr(expr2, GRB.EQUAL, z, "vhStartParkingAtMostone2(" + k + "," + pk + ")");
			}
		} catch (GRBException e) {
			System.out.println("Error code: " + e.getErrorCode() + ". " +
		            e.getMessage());
		}
	}
	
	public void defineEachVehicleArriveDepotAtMostOnceInEachTrip(GRBModel model){
		try{
			for(int k = 0; k < nbVehicles; k++){
				for(int q = 0; q < vhNbTours[k]; q++){
					GRBLinExpr expr1 = new GRBLinExpr();
					GRBLinExpr expr2 = new GRBLinExpr();
					for(ArrayList<Integer> arc : arcVehicles){
						int u = arc.get(0);
						int v = arc.get(1);
						if(centralDepots.contains(u)) {
							String kuvq = k + "-" + u + "-" + v + "-" + q;
							GRBVar x = arc2X.get(kuvq);
							expr1.addTerm(1, x);
						}
						//thoi gian den diem ngay sau central depot
						else if(centralDepots.contains(v)){
							String kuvq = k + "-" + u + "-" + v + "-" + q;
							GRBVar x = arc2X.get(kuvq);
							expr2.addTerm(1, x);
						}
					}
					model.addConstr(expr1, GRB.LESS_EQUAL, 1, "visitOneDepotInATrip1(" + k + "," + q +")");
					model.addConstr(expr2, GRB.LESS_EQUAL, 1, "visitOneDepotInATrip2(" + k + "," + q +")");
				}
			}
		} catch (GRBException e) {
			System.out.println("Error code: " + e.getErrorCode() + ". " +
		            e.getMessage());
		}
	}
	
	public void defineServingTimeConstraintCommon(GRBModel model){
		try{
			for(ArrayList<Integer> arc : arcVehicles){
				for(int k = 0; k < nbVehicles; k++){
					for(int q = 0; q < vhNbTours[k]; q++){
						int u = arc.get(0);
						int v = arc.get(1);
						String kuq = k + "-" + u + "-" + q;
						GRBVar stkuq = point2servingTimeAtPoint.get(kuq);
						
						String kvq = k + "-" + v + "-" + q;
						GRBVar stkvq = point2servingTimeAtPoint.get(kvq);
						
						String kuvq = k + "-" + u + "-" + v + "-" + q;
						GRBVar x = arc2X.get(kuvq);
						
						GRBLinExpr expr1 = new GRBLinExpr();
						expr1.addTerm(1, stkuq);
						expr1.addTerm(-1, stkvq);
						expr1.addTerm(M1, x);
						double tij = travelTime[u][v];
						model.addConstr(expr1, GRB.LESS_EQUAL, M1 - tij - serviceDuration[u], "stTimeAtPointConstr1(" + k + "," + v + "," + q + ")");
					}
				}
			}
		} catch (GRBException e) {
			System.out.println("Error code: " + e.getErrorCode() + ". " +
		            e.getMessage());
		}
	}

	public void defineServingTimeConstraint(GRBModel model){
		try{
			//thoi gian xuat phat parking chuyen dau tien
			for(int k = 0; k < nbVehicles; k++){
				int pk = vh2parking.get(k);
				GRBLinExpr expr1 = new GRBLinExpr();
				String k00 = k + "-" + pk + "-0";
				GRBVar stk00 = point2servingTimeAtPoint.get(k00);
				expr1.addTerm(1, stk00);
				
				model.addConstr(expr1, GRB.EQUAL, earliestArrivalTime[pk], "stAtFirstPk1(" + k + "," + pk + ",0)");
			}
			
			//thoi gian xuat phat tai depot = thoi gian den + thoi gian load hang
			//thoi gian xuat phat tai diem v bang thoi gian tai u cong thoi gian di chuyen + servingTime
			for(int k = 0; k < nbVehicles; k++){
				for(int q = 0; q < vhNbTours[k]; q++){
					for(ArrayList<Integer> arc : arcVehicles){
						int u = arc.get(0);
						int v = arc.get(1);
						if(centralDepots.contains(v))
							continue;
						//thoi gian den diem ngay sau central depot
						if(centralDepots.contains(u)){
							String kuq = k + "-" + u + "-" + q;
							GRBVar stkuq = point2servingTimeAtPoint.get(kuq);
							
							String kvq = k + "-" + v + "-" + q;
							GRBVar stkvq = point2servingTimeAtPoint.get(kvq);
							
							String kuvq = k + "-" + u + "-" + v + "-" + q;
							GRBVar x = arc2X.get(kuvq);
							
							GRBLinExpr expr1 = new GRBLinExpr();
							expr1.addTerm(1, stkuq);
							expr1.addTerm(-1, stkvq);
							expr1.addTerm(M, x);

							for(int p = 0; p < nbProducts; p++){
								String kpq = k + "-" + p + "-" + q;
								GRBVar l = arc2L.get(kpq);
								expr1.addTerm(weightProducts[p]*durationTimeUnit[u], l);
							}
							
							double tij = travelTime[u][v];
							model.addConstr(expr1, GRB.LESS_EQUAL, M - tij - serviceDuration[u], "stTimeAtPointConstr1(" + k + "," + v + "," + q + ")");
							
							GRBLinExpr expr2 = new GRBLinExpr();
							expr2.addTerm(1, stkuq);
							expr2.addTerm(-1, stkvq);
							expr2.addTerm(-M, x);
							
							for(int p = 0; p < nbProducts; p++){
								String kpq = k + "-" + p + "-" + q;
								GRBVar l = arc2L.get(kpq);
								expr2.addTerm(weightProducts[p]*durationTimeUnit[u], l);
							}
							model.addConstr(expr2, GRB.GREATER_EQUAL, -M - tij - serviceDuration[u], "stTimeAtPointConstr2(" + k + "," + v + "," + q + ")");

							continue;
						}
						String kuq = k + "-" + u + "-" + q;
						GRBVar stkuq = point2servingTimeAtPoint.get(kuq);
						
						String kvq = k + "-" + v + "-" + q;
						GRBVar stkvq = point2servingTimeAtPoint.get(kvq);
						
						String kuvq = k + "-" + u + "-" + v + "-" + q;
						GRBVar x = arc2X.get(kuvq);
						
						GRBLinExpr expr1 = new GRBLinExpr();
						expr1.addTerm(1, stkuq);
						expr1.addTerm(-1, stkvq);
						expr1.addTerm(M, x);
						
						double tij = travelTime[u][v];
						model.addConstr(expr1, GRB.LESS_EQUAL, M - tij - serviceDuration[u], "stTimeAtPointConstr1(" + k + "," + v + "," + q + ")");
						
						GRBLinExpr expr2 = new GRBLinExpr();
						expr2.addTerm(1, stkuq);
						expr2.addTerm(-1, stkvq);
						expr2.addTerm(-M, x);
						model.addConstr(expr2, GRB.GREATER_EQUAL, -M - tij - serviceDuration[u], "stTimeAtPointConstr2(" + k + "," + v + "," + q + ")");
					}
				}
			}
			
			//thoi gian den central depot
			for(int k = 0; k < nbVehicles; k++){
				for(int q = 0; q < vhNbTours[k]; q++){
					if(q == 0){
						//chuyen 0: thoi gian den central depot = thoi gian den parking + thoi gian di chuyen
						for(int dp : centralDepots){
							String kdpq = k + "-" + dp + "-" + q;
							GRBVar stkdpq = point2servingTimeAtPoint.get(kdpq);
							int pk = vh2parking.get(k);
							String kpkq = k + "-" + pk + "-" + q;
							GRBVar stkpkq = point2servingTimeAtPoint.get(kpkq);
							
							String kpkdpq = k + "-" + pk + "-" + dp + "-" + q;
							GRBVar x = arc2X.get(kpkdpq);
							
							GRBLinExpr expr1 = new GRBLinExpr();
							expr1.addTerm(1, stkpkq);
							expr1.addTerm(M, x);
							expr1.addTerm(-1, stkdpq);
							
							double tij = travelTime[pk][dp];
							model.addConstr(expr1, GRB.LESS_EQUAL, M - tij - serviceDuration[pk], "stTimeAtDepotConstr1(" + k + "," + dp + "," + q + ")");
							
							GRBLinExpr expr2 = new GRBLinExpr();
							expr2.addTerm(1, stkpkq);
							expr2.addTerm(-M, x);
							expr2.addTerm(-1, stkdpq);
	
							model.addConstr(expr2, GRB.GREATER_EQUAL, -M - tij - serviceDuration[pk], "stTimeAtDepotConstr2(" + k + "," + dp + "," + q + ")");	
						}
					}
					else{
						for(int u : customerPoints){
							for(int dp : centralDepots){
								String kuq1 = k + "-" + u + "-" + (q-1);
								GRBVar stkuq1 = point2servingTimeAtPoint.get(kuq1);
								
								String kdpq = k + "-" + dp + "-" + q;
								GRBVar stkdpq = point2servingTimeAtPoint.get(kdpq);
								
								String kpkdpq1 = k + "-" + u + "-" + dp + "-" + (q-1);
								GRBVar x = arc2X.get(kpkdpq1);
								
								String kq = k + "-" + q;
								GRBVar z = arc2Z.get(kq);
	
								GRBVar w = point2W.get(kuq1);
								
								GRBLinExpr expr1 = new GRBLinExpr();
								expr1.addTerm(1, stkuq1);
								expr1.addTerm(M, z);
								expr1.addTerm(M, w);
//								expr1.addTerm(M, x);
								expr1.addTerm(-1, stkdpq);
								
								double tij = travelTime[u][dp];
								model.addConstr(expr1, GRB.LESS_EQUAL, 2*M - tij - serviceDuration[u], "stTimeAtDepotConstr1(" + k + "," + u + "," + dp + "," + q + ")");
								
								GRBLinExpr expr2 = new GRBLinExpr();
								expr2.addTerm(1, stkuq1);
								expr2.addTerm(-M, z);
								expr2.addTerm(-M, w);
//								expr1.addTerm(-M, x);
								expr2.addTerm(-1, stkdpq);
		
								model.addConstr(expr2, GRB.GREATER_EQUAL, -2*M - tij - serviceDuration[u], "stTimeAtDepotConstr2(" + k + "," + u + "," + dp + "," + q + ")");	
							}
						}
					}
				}
			}
			
//			//thoi gian ve parking
//			for(int k = 0; k < nbVehicles; k++){
//				for(int q = 1; q < vhNbTours[k]; q++){
//					for(int u : truckPoints){
//						if(centralDepots.contains(u))
//							continue;
//						String kuq1 = k + "-" + u + "-" + (q-1);
//						GRBVar stkuq1 = point2servingTimeAtPoint.get(kuq1);
//						
//						for(int dp : centralDepots){
//							String k0q = k + "-" + dp + "-" + q;
//							GRBVar stk0q = point2servingTimeAtPoint.get(k0q);
//							
//							String ku0q1 = k + "-" + u + "-" + dp + "-" + (q-1);
//							GRBVar x = arc2X.get(ku0q1);
//
//							GRBLinExpr expr1 = new GRBLinExpr();
//							expr1.addTerm(1, stkuq1);
//							expr1.addTerm(M, x);
//							expr1.addTerm(-1, stk0q);
//							
//							int tij = travelTime[u][dp];
//							model.addConstr(expr1, GRB.LESS_EQUAL, M - tij - serviceDuration[u], "stTimeAtDepotConstr1(" + k + "," + u + "," + q + ")");
//							
//							GRBLinExpr expr2 = new GRBLinExpr();
//							expr2.addTerm(1, stkuq1);
//							expr2.addTerm(-M, x);
//							expr2.addTerm(-1, stk0q);
//	
//							model.addConstr(expr2, GRB.GREATER_EQUAL, -M - tij - serviceDuration[u], "stTimeAtDepotConstr2(" + k + "," + u + "," + q + ")");	
//						}
//					}
//				}
//			}
			
			
		} catch (GRBException e) {
			System.out.println("Error code: " + e.getErrorCode() + ". " +
		            e.getMessage());
		}
	}
	
	public void defineServingTimeAtPointsWithTimeWindow(GRBModel model){
		try{	
			for(int k = 0; k < nbVehicles; k++){
				for(int q = 0; q < vhNbTours[k]; q++){
					for(int u : truckPoints){
						String kuq = k + "-" + u + "-" + q;
						GRBVar stkuq = point2servingTimeAtPoint.get(kuq);
						
						GRBLinExpr expr1 = new GRBLinExpr();
						expr1.addTerm(1, stkuq);
						for(ArrayList<Integer> arc : arcVehicles){
							int au = arc.get(0);
							int av = arc.get(1);
							if(au == u){
								String kuvq = k + "-" + au + "-" + av + "-" + q;
								GRBVar x = arc2X.get(kuvq);
								expr1.addTerm(M5, x);
							}
						}
						
						model.addConstr(expr1, GRB.LESS_EQUAL, M5 + latestArrivalTime[u], "ServingTimeLELatest(" + k + "," + u + "," + q + ")");
						
						GRBLinExpr expr2 = new GRBLinExpr();
						expr2.addTerm(1, stkuq);
						for(ArrayList<Integer> arc : arcVehicles){
							int au = arc.get(0);
							int av = arc.get(1);
							if(au == u){
								String kuvq = k + "-" + au + "-" + av + "-" + q;
								GRBVar x = arc2X.get(kuvq);
								expr2.addTerm(-M6, x);
							}
						}
						
						model.addConstr(expr2, GRB.GREATER_EQUAL, -M6 + earliestArrivalTime[u], "ServingTimeLEEarliest(" + k + "," + u + "," + q + ")");
					}
				}
			}
		} catch (GRBException e) {
			System.out.println("Error code: " + e.getErrorCode() + ". " +
		            e.getMessage());
		}
	}
	
	public void defineLoadGoodAtDepotConstraint(GRBModel model){		
		try{
			for(int k = 0; k < nbVehicles; k++){
				for(int p = 0; p < nbProducts; p++){
					for(int q = 0; q < vhNbTours[k]; q++){
						GRBLinExpr expr1 = new GRBLinExpr();
						
						String kpq = k + "-" + p + "-" + q;
						GRBVar l = arc2L.get(kpq);
						
						expr1.addTerm(1, l);
						for(ArrayList<Integer> arc : arcVehicles){
							int au = arc.get(0);
							int av = arc.get(1);
							String kuvq = k + "-" + au + "-" + av + "-" + q;
							GRBVar x = arc2X.get(kuvq);
							expr1.addTerm(-demands[au][p], x);
						}
						model.addConstr(expr1, GRB.EQUAL, 0.0, "LoadGoodAtDepot1(" + k + "," + p + "," + q + ")");
					}
				}
			}
		} catch (GRBException e) {
			System.out.println("Error code: " + e.getErrorCode() + ". " +
		            e.getMessage());
		}
	}
	
	public void defineCapacityAtDepotConstraint(GRBModel model){		
		try{
			for(int k = 0; k < nbVehicles; k++){
				for(int q = 0; q < vhNbTours[k]; q++){
					GRBLinExpr expr1 = new GRBLinExpr();
					GRBLinExpr expr2 = new GRBLinExpr();
					for(int p = 0; p < nbProducts; p++){
						String kpq = k + "-" + p + "-" + q;
						GRBVar l = arc2L.get(kpq);
						expr1.addTerm(weightProducts[p], l);
						expr2.addTerm(weightProducts[p], l);
					}
					String kq = k + "-" + q;
					GRBVar z = arc2Z.get(kq);
					expr1.addTerm(M7, z);
					expr2.addTerm(-M7, z);
					model.addConstr(expr1, GRB.LESS_EQUAL, M7 + vhUpperCapacity[k], "upperCapacityAtDepot1(" + k + "," + q + ")");
					model.addConstr(expr2, GRB.GREATER_EQUAL, -M7 + vhLowerCapacity[k], "lowerCapacityAtDepot2(" + k + "," + q + ")");	
				}
			}
		} catch (GRBException e) {
			System.out.println("Error code: " + e.getErrorCode() + ". " +
		            e.getMessage());
		}
	}
	
	public void defineRelatedXZConstraint(GRBModel model){		
		try{
			for(int k = 0; k < nbVehicles; k++){
				for(int q = 0; q < vhNbTours[k]; q++){
					GRBLinExpr expr = new GRBLinExpr();
					
					String kq = k + "-" + q;
					GRBVar z = arc2Z.get(kq);
					expr.addTerm(1, z);
					
					for(ArrayList<Integer> arc : arcVehicles){
						int au = arc.get(0);
						int av = arc.get(1);
						String kuvq = k + "-" + au + "-" + av + "-" + q;
						GRBVar x = arc2X.get(kuvq);
						expr.addTerm(-1, x);
					}
					
					model.addConstr(expr, GRB.LESS_EQUAL, 0.0, "RelatedXZConstr1(" + k + "," + q + ")");
				}
			}
			
			for(int k = 0; k < nbVehicles; k++){
				for(int q = 0; q < vhNbTours[k]; q++){
					for(ArrayList<Integer> arc : arcVehicles){
						GRBLinExpr expr = new GRBLinExpr();
						
						String kq = k + "-" + q;
						GRBVar z = arc2Z.get(kq);
						expr.addTerm(-1, z);
						int au = arc.get(0);
						int av = arc.get(1);
						String kuvq = k + "-" + au + "-" + av + "-" + q;
						GRBVar x = arc2X.get(kuvq);
						expr.addTerm(1, x);
						
						model.addConstr(expr, GRB.LESS_EQUAL, 0.0, "RelatedXZConstr2(" + k + "," + au + "," + av + "," + q + ")");
					}
				}
			}
			
		} catch (GRBException e) {
			System.out.println("Error code: " + e.getErrorCode() + ". " +
		            e.getMessage());
		}
	}
	
	public void defineTripContinuosConstraint(GRBModel model){		
		try{
			for(int k = 0; k < nbVehicles; k++){
				for(int q = 0; q < vhNbTours[k]-1; q++){
					GRBLinExpr expr = new GRBLinExpr();
					String kq = k + "-" + q;
					GRBVar zq = arc2Z.get(kq);
					expr.addTerm(-1, zq);
					
					int nq = q + 1;
					String knq = k + "-" + nq;
					GRBVar znq = arc2Z.get(knq);
					expr.addTerm(1, znq);
					model.addConstr(expr, GRB.LESS_EQUAL, 0.0, "ContinueTrip(" + k + "," + q + ")");
				}
			}
		} catch (GRBException e) {
			System.out.println("Error code: " + e.getErrorCode() + ". " +
		            e.getMessage());
		}
	}
	
	public void defineNextTripConstraint(GRBModel model){		
		try{
//			for(int k = 0; k < nbVehicles; k++){
//				for(int q = 0; q < vhNbTours[k]; q++){
//					for(int i : customerPoints){
//						GRBLinExpr expr = new GRBLinExpr();
//						for(ArrayList<Integer> arc : arcVehicles){
//							int au = arc.get(0);
//							int av = arc.get(1);
//							if(au == i){
//								String kivq = k + "-" + i + "-" + av + "-" + q;
//								GRBVar x1 = arc2X.get(kivq);
//								expr.addTerm(1, x1);
//							}
//							if(av == i){
//								String kivq = k + "-" + au + "-" + i + "-" + q;
//								GRBVar x2 = arc2X.get(kivq);
//								expr.addTerm(-1, x2);
//							}
//						}
//						model.addConstr(expr, GRB.EQUAL, 0.0, "FlowCustomerCustomer(" + k + "," + i + "," + q + ")");
//					}
//				}
//			}
			
			for(int k = 0; k < nbVehicles; k++){
				for(int q = 0; q < vhNbTours[k]-1; q++){
					for(int dp : centralDepots){
						GRBLinExpr expr = new GRBLinExpr();
						for(int i : customerPoints){
							String kidpq = k + "-" + i + "-" + dp + "-" + q;
							String kidpq1 = k + "-" + dp + "-" + i + "-" + (q+1);
							GRBVar x1 = arc2X.get(kidpq);
							GRBVar x2 = arc2X.get(kidpq1);
							expr.addTerm(1, x1);
							expr.addTerm(-1, x2);
						}
						model.addConstr(expr, GRB.EQUAL, 0.0, "FlowCustomerDepot(" + k + "," + dp + "," + q + ")");
					}
				}
			}
			
		} catch (GRBException e) {
			System.out.println("Error code: " + e.getErrorCode() + ". " +
		            e.getMessage());
		}
	}
	
	public void defineLastPointOfTripConstraint(GRBModel model){		
		try{
			for(int k = 0; k < nbVehicles; k++){
				for(int q = 0; q < vhNbTours[k]; q++){
					for(int i : customerPoints){
						GRBLinExpr expr = new GRBLinExpr();
						String kiq = k + "-" + i + "-" + q;
						GRBVar w = point2W.get(kiq);
						expr.addTerm(1, w);
						
						for(int v : outArcVehicles.get(i)){
							if(customerPoints.contains(v)){
								String s = k + "-" + i + "-" + v + "-" + q;
								if(arc2X.get(s) != null){
									GRBVar x = arc2X.get(s);
									expr.addTerm(1, x);
								}
							}
						}
						
						for(int v : inArcVehicles.get(i)){
							String s = k + "-" + v + "-" + i + "-" + q;
							if(arc2X.get(s) != null){
								GRBVar x = arc2X.get(s);
								expr.addTerm(-1, x);
							}
						}

						model.addConstr(expr, GRB.EQUAL, 0.0, "LastPointTrip(" + k + "," + i + "," + q + ")");
					}
				}
			}
		} catch (GRBException e) {
			System.out.println("Error code: " + e.getErrorCode() + ". " +
		            e.getMessage());
		}
	}
	
	public void defineRestrictVehicleCustomerConstraint(GRBModel model){		
		try{
			for(int k = 0; k < nbVehicles; k++){
				for(int i : truckPoints){
					if(vhRestrictCustomers[k][i] == 0){
						GRBLinExpr expr = new GRBLinExpr();
						for(int q = 0; q < vhNbTours[k]; q++){
							for(ArrayList<Integer> arc : arcVehicles){
								int u = arc.get(0);
								int v = arc.get(1);
								if(v == i){
									String kuvq = k + "-" + u + "-" + v + "-" + q;
									GRBVar x = arc2X.get(kuvq);
									expr.addTerm(1, x);
								}
							}
						}
						model.addConstr(expr, GRB.EQUAL, 0.0, "RestrictVehicleCustomer" + k + "," + i + ")");
					}
				}
			}
		} catch (GRBException e) {
			System.out.println("Error code: " + e.getErrorCode() + ". " +
		            e.getMessage());
		}
	}
	
	public void defineRestrictVehicleProductConstraint(GRBModel model){		
		try{
			for(int k = 0; k < nbVehicles; k++){
				for(int p = 0; p < nbProducts; p++){
					if(vhRestrictProducts[k][p] == 0){
						GRBLinExpr expr = new GRBLinExpr();
						for(int q = 0; q < vhNbTours[k]; q++){
							String kpq = k + "-" + p + "-" + q;
							GRBVar l = arc2L.get(kpq);
							expr.addTerm(1, l);
						}
//						for(int q = 0; q < vhNbTours[k]; q++){
//							for(ArrayList<Integer> arc : arcVehicles){
//								int u = arc.get(0);
//								int v = arc.get(1);
//								if(demands[u][p] > 0){
//									String kuvq = k + "-" + u + "-" + v + "-" + q;
//									GRBVar x = arc2X.get(kuvq);
//									expr.addTerm(1, x);
//								}
//							}
//						}
						model.addConstr(expr, GRB.EQUAL, 0.0, "RestrictProductCustomer" + k + "," + p + ")");
					}
				}
			}
		} catch (GRBException e) {
			System.out.println("Error code: " + e.getErrorCode() + ". " +
		            e.getMessage());
		}
	}
	
	public void defineRemainCustomerConstraint(GRBModel model){		
		try{
			for(int k = 0; k < nbVehicles; k++){
				for(int i : customerPoints){
					if(vhRemainCustomers[k][i] == 1){
						GRBLinExpr expr = new GRBLinExpr();
						for(int q = 0; q < vhNbTours[k]; q++){
							for(ArrayList<Integer> arc : arcVehicles){
								int u = arc.get(0);
								int v = arc.get(1);
								if(v == i){
									String kuvq = k + "-" + u + "-" + v + "-" + q;
									GRBVar x = arc2X.get(kuvq);
									expr.addTerm(1, x);
								}
							}
						}
						model.addConstr(expr, GRB.EQUAL, 1.0, "vhRemainCustomer" + k + "," + i + ")");
					}
				}
			}
		} catch (GRBException e) {
			System.out.println("Error code: " + e.getErrorCode() + ". " +
		            e.getMessage());
		}
	}
	
	public void defineServingTimeConstraintsNonLinear(GRBModel model){		
		try{
			//thoi gian xuat phat parking chuyen dau tien
			for(int k = 0; k < nbVehicles; k++){
				int pk = vh2parking.get(k);
				GRBLinExpr expr1 = new GRBLinExpr();
				String k00 = k + "-" + pk + "-0";
				GRBVar stk00 = point2servingTimeAtPoint.get(k00);
				expr1.addTerm(1, stk00);
				
				model.addConstr(expr1, GRB.EQUAL, earliestArrivalTime[pk], "stAtFirstPk1(" + k + "," + pk + ",0)");
			}
			
			
		} catch (GRBException e) {
			System.out.println("Error code: " + e.getErrorCode() + ". " +
		            e.getMessage());
		}
	}
	
	
	public void defineConstraints(GRBModel model){
		flowBalanceConstraint(model);//2
//		flowBalanceConstraintMiller(model);//2
//		flowBalanceOnEachTripConstraint(model);
//		SubtourEliminationConstraint(model);
		SubtourEliminationInTripConstraint(model);
		
		//cong thuc tinh st tai moi diem
		defineServingTimeConstraint(model);//4, 5, 6,7
//		defineServingTimeConstraintCommon(model);
	
		//rang buoc time window cua st tai moi diem
		defineServingTimeAtPointsWithTimeWindow(model);//8,9
//		
		//rang buoc tong so hang load tai depot = tong so hang dem giao
		defineLoadGoodAtDepotConstraint(model);//10
//		
//		//rang buoc can duoi va can tren tai trong cua xe k trong chuyen q
		defineCapacityAtDepotConstraint(model);//11, 12
//		
		//rang buoc moi KH i chi duoc tham 1 lan trong moi xe k moi chuyen q.
		defineEachCustomerVisitedAtMostOnceConstraint(model);//13
//		
		//rang buoc moi xe k chi tham 1 depot trong 1 chuyen q
//		defineEachVehicleArriveDepotAtMostOnceInEachTrip(model);
		
		//rang buoc moi xe k chi xuat phat tu parking 1 lan va quay ve parking 1 lan ///bo
		defineEachVehicleStartsAndComesBackParkingAtMostOnce(model);//14
//		
		//rang buoc moi xe k khong xuat phat tu parking o cac chuyen tiep theo
//		defineEachVehicleNotStartAtParkingForTheNextTrip(model);
//		
		//rang buoc X, Z
		defineRelatedXZConstraint(model);//17,18
//		
//		rang buoc xe k khong duoc xuat phat tai parking khac pk
//		defineVehicleNotStartEndAtOtherParkingConstraint(model);//15,16
//		
		//rang buoc xe k di chuyen q+1 thi phai thuc hien chuyen q
		defineTripContinuosConstraint(model);//19
		
		//rang buoc xe di den dp trong chuyen q thi xe di ra khoi dp trong chuyen q+1
		defineNextTripConstraint(model);//3
		
		//rang buoc tim diem cuoi cung trong chuyen q cua  xe k
		defineLastPointOfTripConstraint(model);
		
		//rang buoc mot so xe ko duoc di den khach i
		defineRestrictVehicleCustomerConstraint(model);
		
		//rang buoc mot so xe ko duoc cho san pham p
		defineRestrictVehicleProductConstraint(model);
		
		//rang buoc xe k phai phuc vu not mot so KH i
		defineRemainCustomerConstraint(model);
	}
	
	public void defineObjective(GRBModel model){
		try{
			GRBLinExpr expr = new GRBLinExpr();
			
			for(ArrayList<Integer> arc : arcVehicles){
				int au = arc.get(0);
				int av = arc.get(1);
				for(int k = 0; k < nbVehicles; k++){
					for(int q = 0; q < vhNbTours[k]; q++){
						String s = k + "-" + au + "-" + av + "-" + q; 
						GRBVar x = arc2X.get(s);
						double d = travelTime[au][av];
						expr.addTerm(d, x);
					}
				}
			}
			model.setObjective(expr, GRB.MINIMIZE);
		} catch (GRBException e) {
			System.out.println("Error code: " + e.getErrorCode() + ". " +
		            e.getMessage());
		}
	}
	
	public HashSet<ArrayList<Integer>> getArcInduced(HashSet<Integer> S, HashSet<ArrayList<Integer>> arcVehicles){
		HashSet<ArrayList<Integer>> A = new HashSet<ArrayList<Integer>>();
		for(ArrayList<Integer> arc : arcVehicles){
			int au = arc.get(0);
			int av = arc.get(1);
			if(S.contains(au) && S.contains(av)){
				ArrayList<Integer> ar = new ArrayList<Integer>();
				ar.add(au);
				ar.add(av);
				A.add(ar);
			}
		}
		return A;
	}
	
	public boolean checkExist(HashSet<ArrayList<Integer>> A, ArrayList<Integer> arc){
		for(ArrayList<Integer> e : A){
			if(e.get(0) == arc.get(0) && e.get(1) == arc.get(1))
				return true;
		}
		return false;
	}
	
	public void calculateServeDurationAtCustomer(){
		serviceDuration = new double[nbCustomers+nbCentralDepots+nbParkings];
		for(int i = 0; i < nbCustomers + nbCentralDepots+nbParkings; i++)
			serviceDuration[i] = 0;
		for(int i = 0; i < nbCustomers + nbCentralDepots+nbParkings; i++){
			serviceDuration[i] += waitingDuration[i];
			for(int p = 0; p < nbProducts; p++)
				serviceDuration[i] += durationTimeUnit[i]*demands[i][p]*weightProducts[p];
		}
	}
	
	public String getNextPoint(int x, int k, int q){
		try{
			for(int i = 0; i < X.size(); i++){
				if(X.get(i).get(GRB.DoubleAttr.X) == 1){
					String[] s = X.get(i).get(GRB.StringAttr.VarName).split(",");
					int kk = Integer.parseInt(s[0].substring(2, s[0].length()));
					int st = Integer.parseInt(s[1]);
					int qq = Integer.parseInt(s[3].substring(0, s[3].length()-1));
					if(st == x && q == qq && k == kk)
						return s[2];
				}
			}
			return null;
			
		} catch (GRBException e) {
			System.out.println("Error code: " + e.getErrorCode() + ". " +
		            e.getMessage());
			return null;
		}
	}
	
	public int getCentralPoint(int k, int q){
		try{
			for(int i = 0; i < X.size(); i++){
				if(X.get(i).get(GRB.DoubleAttr.X) == 1){
					String[] s = X.get(i).get(GRB.StringAttr.VarName).split(",");
					int kk = Integer.parseInt(s[0].substring(2, s[0].length()));
					int st = Integer.parseInt(s[1]);
					int qq = Integer.parseInt(s[3].substring(0, s[3].length()-1));
					if(centralDepots.contains(st) && q == qq && k == kk)
						return st;
				}
			}
			return -1;
			
		} catch (GRBException e) {
			System.out.println("Error code: " + e.getErrorCode() + ". " +
		            e.getMessage());
			return -1;
		}
	}
	
	public void getResult(GRBModel model, String outputFile){
		try{
			double t = System.currentTimeMillis();
			model.optimize();
			
			System.out.println("Optimize done!");
			
			//System.out.println("Obj: " + model.get(GRB.DoubleAttr.ObjVal));
			
			model.update();
			
			for(int i = 0; i < X.size(); i++){
				if(X.get(i).get(GRB.DoubleAttr.X) == 1)
					System.out.println(X.get(i).get(GRB.StringAttr.VarName) + " "
							+ X.get(i).get(GRB.DoubleAttr.X));
			}
			
//			for(int i = 0; i < Y.size(); i++){
//				if(Y.get(i).get(GRB.DoubleAttr.X) == 1)
//					System.out.println(Y.get(i).get(GRB.StringAttr.VarName) + " "
//							+ Y.get(i).get(GRB.DoubleAttr.X));
//			}
			
			for(int i = 0; i < Z.size(); i++){
				if(Z.get(i).get(GRB.DoubleAttr.X) == 1)
					System.out.println(Z.get(i).get(GRB.StringAttr.VarName) + " "
							+ Z.get(i).get(GRB.DoubleAttr.X));
			}
			
			for(int i = 0; i < L.size(); i++){
				//if(L.get(i).get(GRB.DoubleAttr.X) == 1)
					System.out.println(L.get(i).get(GRB.StringAttr.VarName) + " "
							+ L.get(i).get(GRB.DoubleAttr.X));
			}
			
			for(int i = 0; i < W.size(); i++){
				//if(L.get(i).get(GRB.DoubleAttr.X) == 1)
					System.out.println(W.get(i).get(GRB.StringAttr.VarName) + " "
							+ W.get(i).get(GRB.DoubleAttr.X));
			}
			
			
			for(int i = 0; i < servingTimeAtPoint.size(); i++){
				System.out.println(servingTimeAtPoint.get(i).get(GRB.StringAttr.VarName) + " "
						+ servingTimeAtPoint.get(i).get(GRB.DoubleAttr.X));
			}
			
			PrintWriter fw = new PrintWriter(new File(outputFile));
			
			String d = "";
			for(int k = 0; k < nbVehicles; k++){
				for(int q = 0; q < vhNbTours[k]; q++){
					int s = -1;
					HashSet<Integer> visitedNodes = new HashSet<Integer>();
					if(q == 0)
						s = vh2parking.get(k);
					else {
						s = getCentralPoint(k, q);
						visitedNodes.add(s);
					}
					String str = "route[" + k + "]-trip[" + q + "] = " + s + " -> ";
					
					while(true){
						String nextS = getNextPoint(s, k, q);
						if(nextS == null)
							break;
						d += "p1 = " + s + ", p2 = " + nextS 
								+ ", cost = " + travelTime[s][Integer.parseInt(nextS)] + "\n";	
						s = Integer.parseInt(nextS);
						if(Integer.parseInt(nextS) == vh2parking.get(k)
							|| visitedNodes.contains(s)){
							str += nextS;
							break;
						}
						else
							str += nextS + " -> ";
						
						visitedNodes.add(s);
					}
					System.out.println(str);
					fw.println(str);
				}
			}
			
			double runTime = (System.currentTimeMillis() - t)/1000;
			System.out.println("Obj: " + model.get(GRB.DoubleAttr.ObjVal)
					+ ", run time = " + runTime);
			fw.println("Obj: " + model.get(GRB.DoubleAttr.ObjVal)
					+ ", run time = " + runTime);
			System.out.println(d);
			fw.println(d);
			fw.close();
			
		} catch (GRBException e) {
			System.out.println("Error code: " + e.getErrorCode() + ". " +
		            e.getMessage());
		}
		catch(Exception e){
			System.out.println(e);
		}
	}
	
	public static void main(String[] args){
		String dir = "data/MTDLC-VR/";
		String dataFile = dir + "E-1-2parkings-4depots-6vehicles-5customers-MILP.txt";
		String outputFile = dataFile + "-result.txt";
		MILPsolver m = new MILPsolver();

		
		try{
			GRBEnv env   = new GRBEnv();
			env.set("logFile", "vinamilk.log");
			env.start();
			GRBModel model = new GRBModel(env);
			
			m.readData(dataFile);
			
			m.defineVariables(model);
			
			m.defineConstraints(model);
			
			m.defineObjective(model);
			
			System.out.println("Define done!");
			
			m.getResult(model, outputFile);
			
			model.dispose();
			env.dispose();
			
		} catch (GRBException e) {
			System.out.println("Error code: " + e.getErrorCode() + ". " +
		            e.getMessage());
		}
		
	}
}
