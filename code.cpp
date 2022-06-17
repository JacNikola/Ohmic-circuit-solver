#include<bits/stdc++.h>
using namespace std;

//Returns 1 if the node(s) is/are invalid 
bool nodeExceptions(int nodes, int node_1, int node_2) {
	try {
		if(node_1 < 1 || node_1 > nodes) {
			throw node_1;	
		}
	}
	catch(int node_1) {
		cout << "Invalid Node_1\nTry again!\n";	
		return 1;
	}

	try {
		if(node_2 < 1 || node_2 > nodes) {
			throw node_2;	
		}
	}
	catch(int node_2) {
		cout << "Invalid Node_2\nTry again!\n";	
		return 1;
	}

	try {
		if(node_1 == node_2) {
			throw 1;	
		}
	}
	catch(int x) {
		cout << "Node_1 and Node_2 can't be same i.e self loops are not allowed\nTry again!\n";	
		return 1;
	}	

	return 0;
}

//Returns 1 if the circuit is open
bool isCktOpen(vector<vector<int>> &matrix, int n) {
	for(int i = 1; i <= n; i++) {
		int degree = 0;
		for(int j = 1; j <= n; j++) {
			degree += matrix[i][j];
		}
		if(degree < 2) {
			return 1;
		}
	}
	return 0;
} 

//Returns an isolated fundamental set of cycle from a adjacency matrix
bool isolateLoop(int node, int prev, vector<vector<int>> &matrix, vector<int> &ret, vector<bool> &vis, stack<int> &node_stack) {
	vis[node] = 1;
	node_stack.push(node);
	for(int i = 1; i < (int)matrix[0].size(); i++) {
		if(i == prev) continue;
		if(matrix[node][i]) {
			if(vis[i]) {
				ret.push_back(i);
				while(node_stack.top() != i) {
					ret.push_back(node_stack.top());
					node_stack.pop();					
				}
				ret.push_back(node_stack.top());
				return 1;
			}
			else {
				if(isolateLoop(i, node, matrix, ret, vis, node_stack)) return 1;
			}
		}
	}
	if(!node_stack.empty()) node_stack.pop(); 
	return 0;
}

const double EPS = 1e-9;
const int INF = 2; 

int gauss(vector<vector<double>> a, vector<double> & ans) {
    int n = (int) a.size();
    int m = (int) a[0].size() - 1;

    vector<int> where (m, -1);
    for (int col=0, row=0; col<m && row<n; ++col) {
        int sel = row;
        for (int i=row; i<n; ++i)
            if (abs (a[i][col]) > abs (a[sel][col]))
                sel = i;
        if (abs (a[sel][col]) < EPS)
            continue;
        for (int i=col; i<=m; ++i)
            swap (a[sel][i], a[row][i]);
        where[col] = row;

        for (int i=0; i<n; ++i)
            if (i != row) {
                double c = a[i][col] / a[row][col];
                for (int j=col; j<=m; ++j)
                    a[i][j] -= a[row][j] * c;
            }
        ++row;
    }

    ans.assign (m, 0);
    for (int i=0; i<m; ++i)
        if (where[i] != -1)
            ans[i] = a[where[i]][m] / a[where[i]][i];
    for (int i=0; i<n; ++i) {
        double sum = 0;
        for (int j=0; j<m; ++j)
            sum += ans[j] * a[i][j];
        if (abs (sum - a[i][m]) > EPS)
            return 0;
    }

    for (int i=0; i<m; ++i)
        if (where[i] == -1)
            return INF;
    return 1;
}

int main() {

	int nodes;
	cout << "ENTER THE NUMBER OF NODES : ";
	//Each element is sandwiched between two independent nodes
	cin >> nodes;
	if(nodes < 4) {
		cout << "Check carefully if the nodes are numbered properly.\n";
		return 0;
	}

	vector<vector<int>> adj_matrix(nodes+1, vector<int> (nodes+1, 0));
	vector<vector<double>> res_matrix(nodes+1, vector<double> (nodes+1, 0));
	vector<vector<double>> emf_matrix(nodes+1, vector<double> (nodes+1, 0));


	int branches;
	cout << "\nENTER THE NUMBER OF BRANCHES : ";
	cin >> branches;

	//Checks if the circuit is open or closed. Throws ERROR if the circuit is open. 
	try {
		if(branches < nodes) {
			throw branches;	
		}
	}
	catch(int branches) {
		cout << "The circuit is open.";	
		return 0;
	}

	//Branch information
	cout << "***ENTER THE BRANCHES DATA***";
	for(int i = 0; i < branches; i++) {
		int node_1, node_2; 
		cout << "\nBranch " << i+1 << " :\n";
		cout << "Node_1 = ";
		cin >> node_1;
		cout << "Node_2 = ";
		cin >> node_2;
		if(nodeExceptions(nodes, node_1, node_2)) { //If exception caught, i decremented by 1.
			i--;
			continue;
		}
		
		adj_matrix[node_1][node_2] = 1;
		adj_matrix[node_2][node_1] = 1;		
	}	

	if(isCktOpen(adj_matrix, nodes)) { //Checks if there is any node having degree < 2
		cout << "\nERROR: Circuit is open\n";
		return 0;
	}

	vector<vector<int>> adj_matrix_2 = adj_matrix;

	//Resistance information
	int resistors;
	cout << "\nENTER THE NUMBER OF RESISTORS : ";
	cin >> resistors;

	if(resistors) cout << "***ENTER THE RESISTORS DATA***";
	for(int i = 1; i <= resistors; i++) {
		int node_1, node_2;
		double resistance; 
		cout << "\nResistor " << i << " :\n";
		cout << "Node_1 = ";
		cin >> node_1;
		cout << "Node_2 = ";
		cin >> node_2;
		if(nodeExceptions(nodes, node_1, node_2)) { //If exception caught, i decremented by 1.
			i--;
			continue;
		}
		try {
			if(!adj_matrix[node_1][node_2]) {
				throw node_1;
			}
		}
		catch(int node_1) {
			i--;
			cout << "The given branch doesn't exist\nTry Again\n";
			continue;
		}
		cout << "Resistance(ohms) = ";
		cin >> resistance;		

		res_matrix[node_1][node_2] = resistance;
		res_matrix[node_2][node_1] = resistance;		
	}	

	//EMF information
	int emfs; // why double?
	cout << "\nENTER THE NUMBER OF EMFS : ";
	cin >> emfs;

	if(emfs) cout << "***ENTER THE EMFS DATA***";
	for(int i = 1; i <= emfs; i++) {
		int node_1, node_2;
		double emf_val; 
		cout << "\nEMF " << i << " :\n";
		cout << "Node_1(Higher Polarity terminal) = ";
		cin >> node_1;
		cout << "Node_2(Lower Polarity terminal) = ";
		cin >> node_2;
		if(nodeExceptions(nodes, node_1, node_2)) { //If exception caught, i decremented by 1.
			i--;
			continue;
		}
		try {
			if(!adj_matrix[node_1][node_2]) {
				throw node_1;
			}
		}
		catch(int node_1) {
			i--;
			cout << "The given branch doesn't exist\nTry Again\n";
			continue;
		}
		cout << "EMF value (volts) = ";
		cin >> emf_val;		

		emf_matrix[node_1][node_2] = emf_val;
		emf_matrix[node_2][node_1] = -1*emf_val;		
	}	

	//Building spanning tree out of the given graph
	vector<int> vis(nodes+1, 0);
	vector<vector<int>> spanning_matrix(nodes+1, vector<int> (nodes+1, 0));
	vis[1] = 1;

	for(int i = 1; i <= nodes; i++) {
		for(int j = 1; j <= nodes; j++) {
			if(adj_matrix[i][j] == 1 && vis[i] == 1 && vis[j] == 0) {
				spanning_matrix[i][j] = 1;
				spanning_matrix[j][i] = 1;
				vis[j] = 1;				
			}
			if(adj_matrix[i][j] == 1 && vis[i] == 0 && vis[j] == 1) {
				spanning_matrix[i][j] = 1;
				spanning_matrix[j][i] = 1;
				vis[i] = 1;				
			}
		}	
	}

	//Finding fundamental set of cycles
	vector<vector<vector<int>>> loops_vector;

	for(int i = 1; i <= nodes; i++) {
		for(int j = 1; j <= nodes; j++) {
			if(adj_matrix[i][j] ^ spanning_matrix[i][j]) {
				spanning_matrix[i][j] = 1;
				spanning_matrix[j][i] = 1;				
				loops_vector.push_back(spanning_matrix);
				spanning_matrix[i][j] = 0;
				spanning_matrix[j][i] = 0;
				adj_matrix[i][j] = 0;
				adj_matrix[j][i] = 0;
			}	
		}
	}  

	int count_loops = (int)loops_vector.size();
	vector<vector<int>> isolated_loops;

	//Isolating the cycles
	for(auto loop_matrix : loops_vector) {

		vector<bool> vis(nodes+1, 0);
		vector<int> ret;	
		stack<int> node_stack;	
		isolateLoop(1, 0, loop_matrix, ret, vis, node_stack); //marks all the loop nodes as 2 in the adjacency matrix
		for(int i = 1; i <= nodes; i++) {
			for(int j = 1; j <= nodes; j++) {
				loop_matrix[i][j] = 0;
			}
		}

		isolated_loops.push_back(ret);
	}


	//Setting up circuit equations
	map<pair<int, int>, vector<int>> aux_map;
	vector<vector<double>> matrix(count_loops, vector<double> (count_loops+1, 0));

	int loop_num = 0;
	for(auto loop_string : isolated_loops) {  
		for(int i = 0; i < (int)loop_string.size()-1; i++) {
			aux_map[{loop_string[i], loop_string[i+1]}].push_back(loop_num+1);
			aux_map[{loop_string[i+1], loop_string[i]}].push_back(-1 * (loop_num+1));				
			matrix[loop_num][count_loops] += emf_matrix[loop_string[i]][loop_string[i+1]];
		}

		loop_num++;
	}
		
	loop_num = 0;
	for(auto loop_string : isolated_loops) {  
		for(int i = 0; i < (int)loop_string.size()-1; i++) {
			for(auto x : aux_map[{loop_string[i], loop_string[i+1]}]) {
				matrix[loop_num][abs(x)-1] += (x/(abs(x))) * res_matrix[loop_string[i]][loop_string[i+1]];				
			}
		}
		for(auto x : aux_map[{loop_string[(int)loop_string.size()-1], loop_string[0]}]) {
			matrix[loop_num][abs(x)-1] += (x/(abs(x))) * res_matrix[loop_string[(int)loop_string.size()-1]][loop_string[0]];				
		}

		cout << '\n';
		loop_num++;
	}

	vector<double> loop_current;
	gauss(matrix, loop_current);

	cout << "\nCURRENT (in amps)\n\n";
	for(int i = 1; i <= nodes; i++) {
		for(int j = 1; j <= nodes; j++) {
			if(adj_matrix_2[i][j]) {
				cout << i << "->" << j << " = ";
				double c = 0;
				for(auto x : aux_map[{i,j}]) {
					c += ((x/(abs(x)))) * loop_current[abs(x)-1];
				}
				cout << (-1*c) << '\n';
				adj_matrix_2[j][i] = 0;
			}
		}
	}

}
