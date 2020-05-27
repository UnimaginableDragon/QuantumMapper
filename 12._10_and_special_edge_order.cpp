#include <iostream>
#include <queue>
#include <algorithm>
#include <string.h>
#include <set>
#include <climits>
#include <fstream>
#include <bits/stdc++.h>
#define MAX_QUBITS 100
#define INF 9999999999
#define eps 0.000000000000000000001
using namespace std;

struct edge {
	int v1;
	int v2;
	double cost ;
	edge()
	{

	}
	edge(int v1,int v2,double cost)
	{
	    this->v1=v1 ;
	    this->v2=v2 ;
	    this->cost= cost ;
	}
};

struct gate {
	int target;
	int control;
	char type[4];
};

inline bool operator<(const edge& lhs, const edge& rhs) {
	if (lhs.v1 != rhs.v1) {
		return lhs.v1 < rhs.v1;
	}
	return lhs.v2 < rhs.v2;
}

bool edge_cmp(const edge&lhs,const edge&rhs)
{

    if(lhs.cost<rhs.cost+eps){
        return true;
    }
    return false ;


};
struct node {
	double cost_fixed;
	double cost_heur;
	double cost_heur2;
	int depth;
	int* qubits; // get qubit of location -> -1 indicates that there is "no" qubit at a certain location
	int* locations; // get location of qubits -> -1 indicates that a qubit does not have a location -> shall only occur for i > nqubits
	int nswaps;
	int done;
	vector<vector<edge> > swaps;
};

struct node_cmp {
	bool operator()(node& x, node& y) const {
		if ((x.cost_fixed + x.cost_heur + x.cost_heur2) != (y.cost_fixed + y.cost_heur + y.cost_heur2)) {
			return (x.cost_fixed + x.cost_heur + x.cost_heur2) > (y.cost_fixed + y.cost_heur + y.cost_heur2);
		}

		if(x.done == 1) {
			return false;
		}
		if(y.done == 1) {
			return true;
		}

		return x.cost_heur + x.cost_heur2 > y.cost_heur + y.cost_heur2;
	}
};

double** dist;
double** cnot_dist ;
int ** path ;
int positions;
unsigned long ngates = 0;
unsigned int nqubits = 0;


vector <edge> vect_graph ;
std::set<edge> graph;
vector<vector<gate> > layers;
priority_queue<node, std::vector<node>, node_cmp> nodes;
class Qubit{
    int q_num ;
    double gate_error ;
    double readout_errro;
    double T1 , T2 ;
};

class Device{
public:
    int num_qubits ;
    int num_edges ;

    vector < double > qubit_params[MAX_QUBITS] ; //single_gate_error,read_out_error,T1,T2
    vector <pair<int,double>  > adjacency_list[MAX_QUBITS] ;
    std::set<edge> graph ;
    int adj_matrix[MAX_QUBITS][MAX_QUBITS] ;
    Device(){}
    Device(std::ifstream& infile)
    {
        for(int i=0;i<MAX_QUBITS;i++)
        {
            for(int j=0;j<MAX_QUBITS;j++)
                adj_matrix[i][j]=0 ;
        }
        infile >> num_qubits ;
        infile >> num_edges ;
        int q_num ;
        double single_error,readout_error,t1,t2 ;
        for(int i=0;i<num_qubits;i++){
            infile >> q_num >> single_error >> readout_error >> t1 >> t2 ;
            qubit_params[q_num].push_back(single_error);
            qubit_params[q_num].push_back(readout_error);
            qubit_params[q_num].push_back(t1);
            qubit_params[q_num].push_back(t2);
        }




        int v1,v2 ;
        double cx_error ;

        for (int i=0;i<num_edges ;i++)
        {
            infile >>  v1 >> v2 >> cx_error ;
            adjacency_list[v1].push_back(make_pair(v2,cx_error));
            edge e(v1,v2,cx_error) ;
            graph.insert(e);
            vect_graph.push_back(e);

            adj_matrix[v1][v2] = adj_matrix[v2][v1] = 1 ;



        }
       // sort(vect_graph.begin(),vect_graph.end(),edge_cmp);
        //for(int i=0;i<vect_graph.size();i++)
       // {
       //     cout << vect_graph[i].v1 << " " << vect_graph[i].v2 << " " << vect_graph[i].cost << endl;
      //  }

    }
    void print()
    {
        cout << "num_qubits: " << num_qubits << endl;
        for (int i=0;i<num_qubits ;i++)
        {
            cout << qubit_params[i][0] << " " << qubit_params[i][1] << " " << endl;
        }
    }

};

void build_graph(Device & device)
{
    graph.clear();
    set<edge> :: iterator it ;
    positions = device.num_qubits ;
    for (it=device.graph.begin();it!=device.graph.end();it++)
    {
        edge e(it->v1,it->v2,it->cost);
        cout << it->v1 << " " << it->v2 << " " << it->cost << endl;
        graph.insert(e);
    }
}
void build_graph_QX5() {
	graph.clear();
	positions = 16;
	edge e;
	e.v1 = 1;
	e.v2 = 0;
	graph.insert(e);
	e.v1 = 1;
	e.v2 = 2;
	graph.insert(e);
	e.v1 = 2;
	e.v2 = 3;
	graph.insert(e);
	e.v1 = 3;
	e.v2 = 14;
	graph.insert(e);
	e.v1 = 15;
	e.v2 = 0;
	graph.insert(e);
	e.v1 = 15;
	e.v2 = 2;
	graph.insert(e);
	e.v1 = 15;
	e.v2 = 14;
	graph.insert(e);

	e.v1 = 3;
	e.v2 = 4;
	graph.insert(e);

	e.v1 = 5;
	e.v2 = 4;
	graph.insert(e);

	e.v1 = 13;
	e.v2 = 4;
	graph.insert(e);
	e.v1 = 13;
	e.v2 = 14;

	graph.insert(e);
	e.v1 = 12;
	e.v2 = 5;
	graph.insert(e);
	e.v1 = 12;
	e.v2 = 13;

	graph.insert(e);
	e.v1 = 12;
	e.v2 = 11;



	graph.insert(e);
	e.v1 = 6;
	e.v2 = 5;


	graph.insert(e);
	e.v1 = 6;
	e.v2 = 7;



	graph.insert(e);
	e.v1 = 6;
	e.v2 = 11;


	graph.insert(e);
	e.v1 = 11;
	e.v2 = 10;


	graph.insert(e);
	e.v1 = 7;
	e.v2 = 10;

	graph.insert(e);
	e.v1 = 9;
	e.v2 = 10;
	graph.insert(e);

	e.v1 = 9;
	e.v2 = 8;

	graph.insert(e);

	e.v1 = 8;
	e.v2 = 7;

	graph.insert(e);
}


/*bool contains(vector<int,double> v, int e) {
	for (vector<int>::iterator it = v.begin(); it != v.end(); it++) {
		if (it->first == e) {
			return true;
		}
	}
	return false;
}*/


void get_path(int i,int j,vector<int>& arr){
    if (i == j){
        arr.push_back(int(j)) ;
        return ;
    }
    else{
        /// print("i  =",i,"  j=  ",j," ")
        int val = path[int(i)][int(j)] ;
        get_path(i, val, arr) ;
        arr.push_back(int(j)) ;
        return ;
    }
}
double get_max_cnot(vector<int>& p)
{
    double mx =0.0;
    for(int i=0;i<p.size()-1;i++)
    {
        mx=max(mx,cnot_dist[p[i]][p[i+1]]);
    }

    return mx ;
}

void build_dist_table(Device &device) {
	dist = new double*[positions];
	cnot_dist = new double*[positions] ;
    path = new int*[positions] ;

    //cout << " in build dist table 0\n" ;
    for(int i=0;i<positions ;i++)
    {
        dist[i]= new double[positions];
        cnot_dist[i]=new double[positions];
        path[i]= new int[positions];
    }
	for(int i=0;i<positions;i++)
    {
        for(int j=0;j<positions;j++){
            if(i!=j)
                dist[i][j]=INF ;
            else
                dist[i][i]=0 ;

            cnot_dist[i][j] = dist[i][j] ;
        }

    }
    //cout << " in build dist table 1\n" ;
    for (int i=0;i<device.num_qubits;i++)
    {
        for(int j=0;j<device.adjacency_list[i].size();j++)
        {

            int v=device.adjacency_list[i][j].first ;
            double cost=device.adjacency_list[i][j].second ;
            double cost2=cost ;
            if(cost2<cnot_dist[i][v]){
                cnot_dist[i][v]=cost2 ;
            }
            cost2+=2*device.qubit_params[i][0]+2*device.qubit_params[v][0] ;
            if(cost2<cnot_dist[v][i]){
                cnot_dist[v][i]=cost2 ;
            }

            cost = 3 * cost ;
            if(cost<dist[i][v])
            {
                dist[i][v]=cost ;
                path[i][v]=i ;
            }
            if(cost<dist[v][i]){
                dist[v][i]=cost ;
                path[v][i]=v ;
            }
        }
    }

    //cout << " in build dist table 1\n" ;

    //run all pair shortest path algorithm
    //# print(" Now Starting Algorithm")
    for(int k=0;k< positions ;k++){
        for(int i=0;i< positions ; i++){
            for(int  j=0;j< positions ; j++){
                if (dist[i][k] + dist[k][j] < dist[i][j]){
                    dist[i][j] = dist[i][k] + dist[k][j];
                    path[i][j] = path[k][j];
                }
            }
        }
    }



    for(int i=0;i<positions ;i++)
    {
        for(int j=0;j<positions;j++)
        {
            vector<int> p ;
            p.clear();
            get_path(i,j,p);
            dist[i][j]-=get_max_cnot(p);

        }
    }


	for(int i=0;i<vect_graph.size();i++)
    {
        edge e=vect_graph[i];
        double cost1=0;
        double cost2=0;
        for(int i=0;i<positions;i++)
        {
            cost1+=dist[e.v1][i];
            cost2+=dist[e.v2][i];
        }
        e.cost+=cost1+cost2 ;

        vect_graph[i].cost=e.cost ;

    }

    sort(vect_graph.begin(),vect_graph.end(),edge_cmp);
    //for(int i=0;i<p.size();i++)
    //    cout << p[i]<<" -> " ;

   // cout << endl;




}

//A very simplified QASM parser
void read_qasm(std::ifstream& infile) {
	std::string line;

	std::getline(infile, line);
	if(line != "OPENQASM 2.0;") {
		cerr << "ERROR: first line of the file has to be: OPENQASM 2.0;" << endl;
		cerr << line << endl;
		exit(1);
	}

	std::getline(infile, line);
	if(line != "include \"qelib1.inc\";") {
		cerr << "ERROR: second line of the file has to be: include \"qelib1.inc\"" << endl;
		exit(1);
	}

	std::getline(infile, line);
	int n = -1;
	if (sscanf(line.c_str(), "qreg q[%d];", &n) != 1) {
		cerr << "ERROR: failed to parse qasm file: " << line << endl;
		exit(1);
	}
	if (n > positions) {
		cerr << "ERROR: too many qubits for target architecture: " << n << endl;
		exit(2);
	}

	std::getline(infile, line);

	int* last_layer = new int[n];
	for (int i = 0; i < n; i++) {
		last_layer[i] = -1;
	}

	while (std::getline(infile, line)) {

		if (line == "") {
			continue;
		}
		gate g;
		int layer;

		int nq = sscanf(line.c_str(), "%3s q[%d],q[%d];", g.type, &g.control,
				&g.target);

		if (nq == 3) {
			layer = max(last_layer[g.target], last_layer[g.control]) + 1;
			last_layer[g.target] = last_layer[g.control] = layer;
		} else if (nq == 2) {
			g.target = g.control;
			g.control = -1;
			layer = last_layer[g.target] + 1;
			last_layer[g.target] = layer;
		} else {
			double angle;
			if(sscanf(line.c_str(), "rz(%f) q[%d];", &angle, &g.target) == 2) {
				g.control = -1;
				strcpy(g.type, "rz");
				layer = last_layer[g.target] + 1;
				last_layer[g.target] = layer;
			} else {
				cerr << "ERROR: could not read gate: " << line << endl;
				exit(1);
			}
		}
		ngates++;

		if (layers.size() <= layer) {
			layers.push_back(vector<gate>());
		}
		layers[layer].push_back(g);
	}

	nqubits = n;
	delete[] last_layer;
}

void expand_node(const vector<int>& qubits, int qubit, edge *swaps, int nswaps,
		int* used, node base_node, const vector<gate>& gates, double** dist, int next_layer,Device &device) {
   // cout << "in the expand node : next_layer " << next_layer << " qubits " << qubit << endl;
	if (qubit == qubits.size()) {
		//base case: insert node into queue
		//cout << "size is full " << endl;
		if (nswaps == 0) {
			return;
		}
		node new_node;

		new_node.qubits = new int[positions];
		new_node.locations = new int[nqubits];

		memcpy(new_node.qubits, base_node.qubits, sizeof(int) * positions);
		memcpy(new_node.locations, base_node.locations, sizeof(int) * nqubits);
        //cout << "after copying \n "  ;

		new_node.swaps = vector<vector<edge> >();
		new_node.nswaps = base_node.nswaps + nswaps;
		for (vector<vector<edge> >::iterator it2 = base_node.swaps.begin();
				it2 != base_node.swaps.end(); it2++) {
			vector<edge> new_v(*it2);
			new_node.swaps.push_back(new_v);
		}

		new_node.depth = base_node.depth + 5;
		new_node.cost_fixed = base_node.cost_fixed ; /// cost will be calculated from every swap
		new_node.cost_heur = 0;

		vector<edge> new_swaps;
		//cout << " node.cost_fixed: " << new_node.cost_fixed << endl;
		for (int i = 0; i < nswaps; i++) {
			new_swaps.push_back(swaps[i]);
			//cout <<" here I: " << i << " nswaps : " << nswaps << endl;
			int tmp_qubit1 = new_node.qubits[swaps[i].v1];
			int tmp_qubit2 = new_node.qubits[swaps[i].v2];

			new_node.qubits[swaps[i].v1] = tmp_qubit2;
			new_node.qubits[swaps[i].v2] = tmp_qubit1;

			if (tmp_qubit1 != -1) {
				new_node.locations[tmp_qubit1] = swaps[i].v2;
			}
			if (tmp_qubit2 != -1) {
				new_node.locations[tmp_qubit2] = swaps[i].v1;
			}
			//cout <<" before cost update: " << endl;
			//cout <<swaps[i].v1 << " " << swaps[i].v2<< endl;
			new_node.cost_fixed += dist[swaps[i].v1][swaps[i].v2] ; /// dynamic cost
			//cout << "after cost update: " << endl;
		}
		//cout <<" new node cost fixed: " << new_node.cost_fixed << endl;
		new_node.swaps.push_back(new_swaps);
		new_node.done = 1;

		for (vector<gate>::const_iterator it = gates.begin(); it != gates.end();
				it++) {
			gate g = *it;
			if (g.control != -1) {
				new_node.cost_heur = new_node.cost_heur + dist[new_node.locations[g.control]][new_node.locations[g.target]];
                if(new_node.done==1)
                    new_node.done =device.adj_matrix[new_node.locations[g.control]][new_node.locations[g.target]];

			}
			else if(new_node.locations[g.target]!=-1) ///add single qubit gate error.
            {
                new_node.cost_heur+=device.qubit_params[new_node.locations[g.target]][0] ;
            }
		}

		/// adding readout error
		for(int i=0;i<positions;i++)
        {
            if(new_node.qubits[i]!=-1)new_node.cost_heur+=device.qubit_params[i][1] ;
        }

		//Calculate heuristics for the cost of the following layer
		new_node.cost_heur2 = 0;
		if(next_layer != -1) {
			for (vector<gate>::const_iterator it = layers[next_layer].begin(); it != layers[next_layer].end();
							it++) {
				gate g = *it;
				if (g.control != -1) {
					if(new_node.locations[g.control] == -1 && new_node.locations[g.target]) {
					//ignore this case
					} else if(new_node.locations[g.control] == -1) {
						double min = INF;
						for(int i=0; i< positions; i++) {
							if(new_node.qubits[i] == -1 && dist[i][new_node.locations[g.target]] < min) {
								min = dist[i][new_node.locations[g.target]];
							}
						}
						new_node.cost_heur2 = new_node.cost_heur2 + min;
					} else if(new_node.locations[g.target] == -1) {
						double min = INF;
						for(int i=0; i< positions; i++) {
							if(new_node.qubits[i] == -1 && dist[new_node.locations[g.control]][i] < min) {
								min = dist[new_node.locations[g.control]][i];
							}
						}
						new_node.cost_heur2 = new_node.cost_heur2 + min;
					} else {
						new_node.cost_heur2 = new_node.cost_heur2 + dist[new_node.locations[g.control]][new_node.locations[g.target]];
					}
				}
				else if(g.control == -1 && new_node.locations[g.target]!=-1) {
					//ignore this case
					/// add single qubit gate error
					new_node.cost_heur2+=device.qubit_params[new_node.locations[g.target]][0] ;
                }
			}
		}

		nodes.push(new_node);
	} else {
	    //cout << "in expand of next layer : " << qubit << endl;
		expand_node(qubits, qubit + 1, swaps, nswaps, used, base_node, gates,
				dist, next_layer,device);
       // cout <<" after expansion of next layer: " << qubit << endl;

		for (vector<edge>::iterator it = vect_graph.begin(); it != vect_graph.end(); it++) {
			edge e = *it;

			//cout << e.v1 <<"  "<<e.v2 << "  " << e.cost << endl;
			if (e.v1 == base_node.locations[qubits[qubit]]
					|| e.v2 == base_node.locations[qubits[qubit]]) {
				if (!used[e.v1] && !used[e.v2]) {
					used[e.v1] = 1;
					used[e.v2] = 1;
					swaps[nswaps].v1 = e.v1;
					swaps[nswaps].v2 = e.v2;
					//cout << "iner if " << endl;
					expand_node(qubits, qubit + 1, swaps, nswaps + 1, used,
							base_node, gates, dist, next_layer,device);
                    //cout << "after inner if" << endl;
					used[e.v1] = 0;
					used[e.v2] = 0;
				}
			}
		}
	}
}

int getNextLayer(int layer) {
	int next_layer = layer+1;
	while(next_layer < layers.size()) {
		for(vector<gate>::iterator it = layers[next_layer].begin(); it != layers[next_layer].end(); it++) {
			if(it->control != -1) {
				return next_layer;
			}
		}
		next_layer++;
	}
	return -1;
}

node a_star_fixlayer(int layer, int* map, int* loc, double** dist,Device &device ) {

	int next_layer = getNextLayer(layer);
    //cout << " in a_star_fix layer :  " << layer << endl;
	node n;
	n.cost_fixed = 0;
	n.cost_heur = n.cost_heur2 = 0;
	n.qubits = new int[positions];
	n.locations = new int[nqubits];
	n.swaps = vector<vector<edge> >();
	n.done = 1;

	vector<gate> v = vector<gate>(layers[layer]);
	vector<int> considered_qubits;

	//Find a mapping for all logical qubits in the CNOTs of the layer that are not yet mapped

	for (vector<gate>::iterator it = v.begin(); it != v.end(); it++) {
		gate g = *it;
		if (g.control != -1) {
			considered_qubits.push_back(g.control);
			considered_qubits.push_back(g.target);
			if(loc[g.control] == -1 && loc[g.target] == -1) {
				vector<edge> possible_edges;
				for(vector<edge>::iterator it = vect_graph.begin(); it != vect_graph.end(); it++) {
					if(map[it->v1] == -1 && map[it->v2] == -1) {
						possible_edges.push_back(*it);
					}
				}
				if(!possible_edges.empty()) {
					edge e = *possible_edges.begin();
                    cout << e.v1 <<"  "<<e.v2 << "  " << e.cost << endl;

					loc[g.control] = e.v1;
					map[e.v1] = g.control;
					loc[g.target] = e.v2;
					map[e.v2] = g.target;
				} else {
					cout << "no edge available!";
					exit(1);
				}
			} else if(loc[g.control] == -1) {
				double min = INF;
				int min_pos = -1;
				for(int i=0; i< positions; i++) {
					if(map[i] == -1 && dist[i][loc[g.target]] < min) {
						min = dist[i][loc[g.target]];
						min_pos = i;
					}
				}
				map[min_pos] = g.control;
				loc[g.control] = min_pos;
			} else if(loc[g.target] == -1) {
				double min = INF;
				int min_pos = -1;
				for(int i=0; i< positions; i++) {
					if(map[i] == -1 && dist[loc[g.control]][i] < min) {
						min = dist[loc[g.control]][i];
						min_pos = i;
					}
				}
				map[min_pos] = g.target;
				loc[g.target] = min_pos;
			}
			if(n.cost_heur<dist[loc[g.control]][loc[g.target]])
            {
                n.cost_heur = dist[loc[g.control]][loc[g.target]];

            }
            if(n.done==1)
                    n.done=device.adj_matrix[loc[g.control]][loc[g.target]] ;

		} else if(g.control== -1 && loc[g.target]!=-1) {
					//ignore this case
					/// add single qubit gate error
            n.cost_heur2+=device.qubit_params[loc[g.target]][0] ;
        }
	}

	//if(n.cost_heur > 4) {
	//	n.done = 0;
	//}


	memcpy(n.qubits, map, sizeof(int) * positions);
	memcpy(n.locations, loc, sizeof(int) * nqubits);

/// adding readout error
    for(int i=0;i<positions;i++)
    {
        if(n.qubits[i]!=-1)n.cost_fixed+=device.qubit_params[i][1] ;
    }
    nodes.push(n);


	int *used = new int[positions];
	for (int i = 0; i < positions; i++) {
		used[i] = 0;
	}
	edge *edges = new edge[considered_qubits.size()];
    //cout << "after mapping  start  A*\n" ;
	//Perform an A* search to find the cheapest permuation
	while (!nodes.top().done) {
		node n = nodes.top();
		nodes.pop();

		expand_node(considered_qubits, 0, edges, 0, used, n, v, dist, next_layer,device);

		delete[] n.locations;
		delete[] n.qubits;
	}

	node result = nodes.top();
	nodes.pop();

	//clean up
	delete[] used;
	delete[] edges;
	while (!nodes.empty()) {
		node n = nodes.top();
		nodes.pop();
		delete[] n.locations;
		delete[] n.qubits;
	}
	return result;
}

int main(int argc, char** argv) {
    Device device ;
	if(argc != 5) {
		cout << "Usage: ./imb_mapping <device_conf_file> <input_file> <output_file> <output_layout_file>" << endl;
		exit(0);
	}

	//build_graph_QX3();
    std::ifstream infile ;
	if(strcmp(argv[1],"null")==0)
    {
        build_graph_QX5();
        //build_dist_table(graph);

        cout << " building default: " << endl;
    }
    else {
        //cout << "before taking input \n" ;
        infile.open(argv[1]);
        Device dev(infile);
        //cout << "after taking input \n" ;
        dev.print();
        //cout << "after printing device info\n" ;
        infile.close();
        build_graph(dev);
        //cout << "after building device graph\n" ;
        build_dist_table(dev);
        //cout << "after building dist_table \n " ;

        device=dev ;

    }




	infile.open(argv[2]);
	read_qasm(infile);
	infile.close();
 //cout << "2 " << endl ;

	unsigned int width = 0;
	for (vector<vector<gate> >::iterator it = layers.begin(); it != layers.end(); it++) {
		if ((*it).size() > width) {
			width = (*it).size();
		}
	}

	//char* bName = basename(argv[1]);
	char* bName = (argv[2]);
	cout << "Circuit name: " << bName << " (requires " << nqubits << " qubits)" << endl;

	cout << endl << "Before mapping: " << endl;
	cout << "  elementary gates: " << ngates << endl;
	cout << "  depth: " << layers.size() << endl;

	int *locations = new int[nqubits];
	int *qubits = new int[positions];

	for (int i = 0; i < positions; i++) {
		qubits[i] = -1;
	}
	for(int i = 0; i < nqubits; i++) {
		locations[i] = qubits[i] = i;
	}

	//Start mapping algorithm
	clock_t begin_time = clock();
     //cout <<"debug: 3 " << endl;
	//Initially, no physical qubit is occupied
	for (int i = 0; i < positions; i++) {
			qubits[i] = -1;
	}

	//Initially, no logical qubit is mapped to a physical one
	for(int i = 0; i < nqubits; i++) {
		locations[i] = -1;
	}

	vector<gate> all_gates;
	int total_swaps = 0;
    //cout <<"debug: 4 " << endl;
	//Fix the mapping of each layer
	for (int i = 0; i < layers.size(); i++) {
        //cout << "layer no: " << i << endl;
		node result = a_star_fixlayer(i, qubits, locations, dist,device);
        //cout << "after fixing maping  \n" << endl;
		delete[] locations;
		delete[] qubits;
		locations = result.locations;
		qubits = result.qubits;

		vector<gate> h_gates = vector<gate>();
        //cout << "debug: 5 -> " << i << endl;
		//The first layer does not require a permutation of the qubits
		if (i != 0) {
			//Add the required SWAPs to the circuits
			for (vector<vector<edge> >::iterator it = result.swaps.begin();
					it != result.swaps.end(); it++) {
				for (vector<edge>::iterator it2 = it->begin(); it2 != it->end();
						it2++) {

					edge e = *it2;
					gate cnot;
					gate h1;
					gate h2;
					if (graph.find(e) != graph.end()) {
						cnot.control = e.v1;
						cnot.target = e.v2;
					} else {
						cnot.control = e.v2;
						cnot.target = e.v1;

						int tmp = e.v1;
						e.v1 = e.v2;
						e.v2 = tmp;
						if (graph.find(e) == graph.end()) {
							cerr << "ERROR: invalid SWAP gate" << endl;
							exit(2);
						}
					}
					strcpy(cnot.type, "cx");
					strcpy(h1.type, "h");
					strcpy(h2.type, "h");
					h1.control = h2.control = -1;
					h1.target = e.v1;
					h2.target = e.v2;

					gate gg;
					gg.control = cnot.control;
					gg.target = cnot.target;
					strcpy(gg.type, "SWP");

					all_gates.push_back(cnot);
					all_gates.push_back(h1);
					all_gates.push_back(h2);
					all_gates.push_back(cnot);
					all_gates.push_back(h1);
					all_gates.push_back(h2);
					all_gates.push_back(cnot);
					//Insert a dummy SWAP gate to allow for tracking the positions of the logical qubits
					all_gates.push_back(gg);
					total_swaps++;
				}
			}
		}

		//Add all gates of the layer to the circuit
		//cout <<" debug: 6 " << endl;
		vector<gate> layer_vec = layers[i];
		for (vector<gate>::iterator it = layer_vec.begin();
				it != layer_vec.end(); it++) {
                //cout << " debug : 7 " << endl;
			gate g = *it;
			if (g.control == -1) {
				//single qubit gate
				//cout <<" debug : 8 : " << endl;
				if(locations[g.target] == -1) {
                        //cout <<"debug: 9 " << endl;
					//handle the case that the qubit is not yet mapped. This happens if the qubit has not yet occurred in a CNOT gate
					gate g2 = g;
					g2.target = -g.target -1;
					all_gates.push_back(g2);
                //cout <<"debug: 10 " <<endl;
				} else {
					//Add the gate to the circuit
					g.target = locations[g.target];
					all_gates.push_back(g);
				}
			} else {
				//CNOT gate
				g.target = locations[g.target];
				g.control = locations[g.control];

				edge e;
				e.v1 = g.control;
				e.v2 = g.target;

				if (graph.find(e) == graph.end()) {
					//flip the direction of the CNOT by inserting H gates
					e.v1 = g.target;
					e.v2 = g.control;
					if (graph.find(e) == graph.end()) {
						cerr << "ERROR: invalid CNOT: " << e.v1 << " - " << e.v2
								<< endl;
						exit(3);
					}
					gate h;
					h.control = -1;
					strcpy(h.type, "h");
					h.target = g.target;
					all_gates.push_back(h);

					h_gates.push_back(h);
					h.target = g.control;
					all_gates.push_back(h);

					h_gates.push_back(h);
					int tmp = g.target;
					g.target = g.control;
					g.control = tmp;
				}
				all_gates.push_back(g);
			}
			//cout << "debug : 11 " << endl;
		}
		//cout << "debug : 12 " << endl;
        //cout << " h_gates.size(): "  << h_gates.size() << endl ;
		if (h_gates.size() != 0) {
                //cout <<"debug: 13 " << endl;
			if (result.cost_heur == 0) {
				cerr << "ERROR: invalid heuristic cost!" << endl;
				exit(2);
			}

			for (vector<gate>::iterator it = h_gates.begin();
					it != h_gates.end(); it++) {
				all_gates.push_back(*it);
			}
		}
		//cout << "after debug : 12 " << endl;

	}
	int * final_locations= new int[nqubits];
	for(int i=0;i<nqubits;i++)
    {
        final_locations[i]=locations[i];
    }
	//cout <<"debug: 5 " <<endl;

	//Fix the position of the single qubit gates
	for(vector<gate>::reverse_iterator it = all_gates.rbegin(); it != all_gates.rend(); it++) {
		if(strcmp(it->type, "SWP") == 0) {
			int tmp_qubit1 = qubits[it->control];
			int tmp_qubit2 = qubits[it->target];
			qubits[it->control] = tmp_qubit2;
			qubits[it->target] = tmp_qubit1;

			if(tmp_qubit1 != -1) {
				locations[tmp_qubit1] = it->target;
			}
			if(tmp_qubit2 != -1) {
				locations[tmp_qubit2] = it->control;
			}
		}
		if(it->target < 0) {
			int target = -(it->target +1);

			if(locations[target] == -1) {
				//This qubit occurs only in single qubit gates -> it can be mapped to an arbirary physical qubit
				int loc = 0;
				while(qubits[loc] != -1) {
					loc++;
				}
				locations[target] = loc;
				qubits[loc]=target ;
			}
            it->target = locations[target];

		}
	}

    for(int i=0;i<nqubits;i++)
    {
        if(final_locations[i]==-1)
            final_locations[i]=locations[i];
    }
	int *last_layer = new int[positions];
	for(int i=0; i<positions; i++) {
		last_layer[i] = -1;
	}

	vector<vector<gate> > mapped_circuit;


	//build resulting circuit
	for(vector<gate>::iterator it = all_gates.begin(); it != all_gates.end(); it++) {
		if(strcmp(it->type, "SWP") == 0) {
			continue;
		}
		if(it->control == -1) {
			//single qubit gate
			gate g = *it;
			int layer = last_layer[g.target] + 1;

			if (mapped_circuit.size() <= layer) {
				mapped_circuit.push_back(vector<gate>());
			}
			mapped_circuit[layer].push_back(g);
			last_layer[g.target] = layer;
		} else {
			gate g = *it;
			int layer = max(last_layer[g.control], last_layer[g.target]) + 1;
			if (mapped_circuit.size() <= layer) {
				mapped_circuit.push_back(vector<gate>());
			}
			mapped_circuit[layer].push_back(g);

			last_layer[g.target] = layer;
			last_layer[g.control] = layer;
		}
	}

    //cout << "4" << endl;
	double time = double(clock() - begin_time) / CLOCKS_PER_SEC;

	cout << endl << "After mapping (no post mapping optimizations are conducted): " << endl;
	cout << "  elementary gates: " << all_gates.size()-total_swaps << endl;
	cout << "  depth: " << mapped_circuit.size() << endl;

	cout << endl << "The mapping required " << time << " seconds" << endl;

	//cout << endl << "Initial mapping of the logical qubits (q) to the physical qubits (Q) of the IBM QX5 architecture: " << endl;
    ofstream lof(argv[4]);
	for(int i=0; i<nqubits; i++) {
		//cout << "  q" << i << " is initially mapped to Q" << locations[i] << endl;
		lof << i << " " << final_locations[i] << endl;
	}
	lof<<time<<endl;

	//Dump resulting circuit

	ofstream of(argv[3]);

	of << "OPENQASM 2.0;" << endl;
	of << "include \"qelib1.inc\";" << endl;
	of << "qreg q["<<positions<<"];" << endl;
	of << "creg c["<<positions<<"];" << endl;

	for (vector<vector<gate> >::iterator it = mapped_circuit.begin();
			it != mapped_circuit.end(); it++) {
		vector<gate> v = *it;
		for (vector<gate>::iterator it2 = v.begin(); it2 != v.end(); it2++) {
			of << it2->type << " ";
			if (it2->control != -1) {
				of << "q[" << it2->control << "],";
			}
			of << "q[" << it2->target << "];" << endl;
		}
	}

	delete[] final_locations;
	delete[] locations;
	delete[] qubits;
	delete[] last_layer;

	return 0;
}
