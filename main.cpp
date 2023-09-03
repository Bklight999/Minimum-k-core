#include <fstream>
#include <iostream>
#include <bitset>
#include<queue>
#include <vector>
#include <set>
#include <ctime>
#include<cmath>
#include<string.h>
#include <time.h>
#include<algorithm>

#define FI first
#define SE second

using namespace std;



int nnodes, nedges, K, S;//节点数、边数、参数K、S

const int MaxN = 3000000; // maximum number of the vertices
const int MaxM = 20000000; // maximum number of the edges
const int Maxn = 800000; // maximum number of the vertices after the reduction

double timelimit = 8640000;





//std::pair <int, int> E[MaxM], newE[MaxM];;//<起点，终点>
//
//int* DEG, * head, * nxt, * to, * tempdeg;
//int startfromzero, n, ecnt;
//int* que, * nV, * oriID;
//bool* outcore, * tempoutcore;
//int* degree, * degree_temp;
//vector<double> pro(MaxN, 0);
//vector<double> reversedeg(MaxN, -1.0);
//bool* del, * ins;
//int* notadj, * adj;
//std::vector <int> delvex, svex;
//int LB, UB;
//int finaltime;
//int* exist_oldid, * triangles;
//vector<string> fileset;

std::pair <int, int> E[MaxM], newE[MaxM], temp_E[MaxM];//<起点，终点>

int DEG[MaxN], head[MaxN], nxt[2 * MaxM], to[2 * MaxM], ecnt, tempdeg[MaxN];
double vertex_weight[MaxN];
int startfromzero;

int que[MaxN], nV[MaxN], oriID[MaxN];

bool outcore[MaxN];
bool tempoutcore[MaxN];

int degree[MaxN], n, degree_temp[MaxN];
//std::bitset <Maxn> Graph[Maxn];

vector<double> pro(MaxN, 0);
vector<double> reversedeg(MaxN, -1.0);

bool del[MaxN], ins[MaxN];

int notadj[Maxn], adj[Maxn], tempadj[Maxn];

std::vector <int> delvex, svex;


int LB, UB;
int finaltime;

int exist_oldid[MaxN];
int triangles[2 * MaxM];
vector<string> fileset;
int* appear;

//void initialize()
//{
//    DEG = new int[nnodes];
//    head = new int[nnodes];
//    nxt = new int[2 * nedges];
//    to = new int[2 * nedges];
//    tempdeg = new int[nnodes];
//    que = new int[nnodes];
//    nV = new int[nnodes];
//    oriID = new int[nnodes];
//    outcore = new bool[nnodes];
//    tempoutcore = new bool[nnodes];
//    degree = new int[nnodes];
//    degree_temp = new int[nnodes];
//    del = new bool[nnodes];
//    ins = new bool[nnodes];
//    notadj = new int[nnodes];
//    adj = new int[nnodes];
//    exist_oldid = new int[nnodes];
//    triangles = new int[2 * nedges];
//}



vector<std::string> readtxt_name(string file)
{
    ifstream infile;
    infile.open(file.data());
    cout << "infile.is_open() is : " << infile.is_open() << endl;

    if (infile.is_open() == 0) {
        cout << "file not existed" << endl;
    }
    vector< std::string> txt_name;
    string str;

    while (getline(infile, str))
    {
        if (str == "EOF")break;
        txt_name.push_back(str);

    }
    infile.close();
    return txt_name;
}

void refresh()
{
    for (int i = 0; i < MaxN; i++)
    {
        DEG[i] = 0;
        head[i] = 0;
        tempdeg[i] = 0;
        que[i] = 0;
        nV[i] = 0;
        oriID[i] = 0;
        outcore[i] = false;
        tempoutcore[i] = false;
        degree[i] = 0;
        degree_temp[i] = 0;
        pro[i] = 0;
        reversedeg[i] = -1.0;
        del[i] = false;
        ins[i] = false;
    }
    for (int i = 0; i < Maxn; i++)
    {
        notadj[i] = 0;
        adj[i] = 0;
    }
    for (int i = 0; i < MaxM; i++)
    {
        E[i].first = 0;
        E[i].second = 0;
    }
    for (int i = 0; i < 2 * MaxM; i++)
    {
        nxt[i] = 0;
        to[i] = 0;
    }
    n = 0;
    ecnt = 0;
    delvex.clear();
    svex.clear();
    LB = UB = 0;
}

void output_max_k_core()
{
    vector<int> ans_v;
    for (int i = 0; i < nnodes; i++)
        if (!outcore[i]) ans_v.push_back(i);
    ofstream outfile;
    outfile.open("/home/zhangqifan/min_k_core/max_k_core.txt");
    set<pair<int, int>> ans_e;
    for (auto v : ans_v)
    {
        for (int e = head[v]; ~e; e = nxt[e])
        {
            int u = to[e];
            if (!outcore[u])
            {
                if (v < u)
                {
                    ans_e.insert(make_pair(v, u));
                }
                else
                {
                    ans_e.insert(make_pair(u, v));
                }
            }
        }
    }
    for (auto e : ans_e)
        outfile << e.first << " " << e.second << " " << endl;
    outfile.close();
}

bool inans(int u, vector<int> ans)
{
    for (int i = 0; i < ans.size(); i++)
    {
        if (ans[i] == u) return true;
    }
    return false;
}

void output_min_k_core()
{
    vector<int> ans_v;
    int cnt = 0;
    for (int i = 0; i < n; i++)
    {
        if (!del[i] && cnt < LB)
        {
            cnt++;
            ans_v.push_back(i);
        }
    }
    ofstream outfile;
    outfile.open("/home/zhangqifan/min_k_core/min_k_core.txt");
    set<pair<int, int>> ans_e;
    for (auto v : ans_v)
    {
        for (int e = head[v]; ~e; e = nxt[e])
        {
            int u = to[e];
            if (!del[u] && inans(u, ans_v))
            {
                if (v < u)
                {
                    ans_e.insert(make_pair(oriID[v], oriID[u]));
                }
                else
                {
                    ans_e.insert(make_pair(oriID[u], oriID[v]));
                }
            }
        }
    }
    for (auto e : ans_e)
        outfile << e.first << " " << e.second << " " << endl;
    outfile.close();
}

bool check_neighbor(int u, int v)
{
    for (int e = head[u]; ~e; e = nxt[e])
    {
        if (to[e] == v)
        {
            return true;
        }
    }
    return false;
}

void addedge(int* h, int v) {
    nxt[ecnt] = *h;
    to[ecnt] = v;
    *h = ecnt++;
}

bool checkedge(int u, int v) {
    for (int e = head[u]; ~e; e = nxt[e]) {
        if (to[e] == v)
            return true;
    }
    return false;
}

//检测图是否为空
bool checkempty(bool outcore[])
{
    //注意这里是i<n而不是i<exist.size()
    for (int i = 0; i < nnodes; i++)
        if (!outcore[i])
            return false;
    return true;
}

void clean()
{
    queue<int> waitdel;
    for (int i = 0; i < nnodes; i++)
    {
        if (DEG[i] < K && !outcore[i])
        {
            outcore[i] = true;
            DEG[i] = 0;
            waitdel.push(i);
        }
    }
    while (!waitdel.empty())
    {
        int temp = waitdel.front();
        waitdel.pop();
        for (int e = head[temp]; ~e; e = nxt[e])
        {
            int v = to[e];
            if (!outcore[v] && --DEG[v] < K)
            {
                outcore[v] = true;
                DEG[v] = 0;
                waitdel.push(v);
            }
        }
    }
    /*
    int fr = 0, re = 0;
    for (int i = 0; i < nnodes; ++i) {
        if (DEG[i] < K) {
            que[re++] = i;
            outcore[i] = true;//outcore[i]表示是否删除
        }
    }

    while (fr ^ re) {
        int u = que[fr++];
        for (int e = head[u]; ~e; e = nxt[e]) {
            int v = to[e];
            if (outcore[v]) continue;
            if (--DEG[v] < K) {
                que[re++] = v;
                outcore[v] = true;
            }
        }
    }//不断重复删除
    */
}

void oldpreprocess()
{
    memset(DEG, 0, nnodes * sizeof(int));
    memset(head, -1, nnodes * sizeof(int));
    ecnt = 0;
    for (int i = 0; i < nedges; ++i) {
        ++DEG[E[i].FI];//边的两端度数加1
        ++DEG[E[i].SE];
        addedge(head + E[i].FI, E[i].SE);
        addedge(head + E[i].SE, E[i].FI);
    }
    memset(outcore, 0, nnodes * sizeof(bool));
    memset(exist_oldid, 1, sizeof(int) * nnodes);
    double total = 0;
    for (int i = 0; i < nnodes; i++)
        total += DEG[i];
    cout << "average degree is" << total / nnodes << endl;
    int fr = 0, re = 0;
    for (int i = 0; i < nnodes; ++i) {
        if (DEG[i] < K) {
            que[re++] = i;
            outcore[i] = true;//outcore[i]表示是否删除
        }
    }

    while (fr ^ re) {
        int u = que[fr++];
        exist_oldid[u] = 0;
        for (int e = head[u]; ~e; e = nxt[e]) {
            int v = to[e];
            if (outcore[v]) continue;
            if (--DEG[v] < K) {
                que[re++] = v;
                outcore[v] = true;
            }
        }
    }//不断重复删除
    output_max_k_core();
    n = 0;
    for (int i = 0; i < nnodes; ++i) {
        if (!outcore[i])
            nV[i] = n, vertex_weight[nV[i]] = vertex_weight[i], oriID[n++] = i;
    }
    //n为预处理后的节点数量
    //for(int i = 0; i < n; ++ i) Graph[i].reset();
    memset(head, -1, n * sizeof(int));
    memset(degree, 0, n * sizeof(int));
    memset(degree_temp, 0, n * sizeof(int));

    ecnt = 0;
    for (int i = 0; i < nedges; ++i) {
        int u = E[i].FI, v = E[i].SE;
        if (outcore[u] || outcore[v]) continue;
        u = nV[u]; v = nV[v];
        ++degree[u]; ++degree[v];
        ++degree_temp[u]; ++degree_temp[v];
        addedge(head + u, v);
        addedge(head + v, u);
    }

    cout << "预处理后的点数量为" << n << endl;
    memset(outcore, 0, sizeof(bool) * n);

}

void preprocess() {
    //Start to preprocess the graph    
    memset(DEG, 0, nnodes * sizeof(int));
    memset(head, -1, nnodes * sizeof(int));
    ecnt = 0;
    for (int i = 0; i < nedges; ++i) {
        ++DEG[E[i].FI];//边的两端度数加1
        ++DEG[E[i].SE];
        addedge(head + E[i].FI, E[i].SE);
        addedge(head + E[i].SE, E[i].FI);
    }
    memset(outcore, 0, nnodes * sizeof(bool));
    //preprocess
    double c = 0.0, resum = 0.0;

    while (!checkempty(outcore))
    {
        fill(pro.begin(), pro.begin() + nnodes, 0);
        fill(reversedeg.begin(), reversedeg.begin() + nnodes, -1.0);
        c = 0.0;
        resum = 0.0;
        for (int i = 0; i < nnodes; i++)
        {
            tempdeg[i] = DEG[i];
            tempoutcore[i] = outcore[i];
        }
        for (int i = 0; i < nnodes; i++)
        {
            if (DEG[i] > 0 && !outcore[i])
                reversedeg[i] = 1 / (DEG[i] * 1.0);
            //deg[i]等于0该怎么办呢？deg[i]为0一定会被删除，因为1-core也至少度数为1
            else
            {
                reversedeg[i] = -1.0;
                outcore[i] = true;
                DEG[i] = 0;
            }
        }

        for (int i = 0; i < nnodes; i++)
        {
            if (reversedeg[i] != -1.0 && !outcore[i])//等于-1.0说明度数为0，直接删除就好
                resum += reversedeg[i];
        }
        c = 1 / resum;
        //计算每个节点的probability，这里用的是local degree weighting
        for (int i = 0; i < nnodes; i++)
        {
            if (DEG[i] > 0 && !outcore[i])
                pro[i] = c / (DEG[i] * 1.0);
            //避免重复删除已经删除的节点
            else
                pro[i] = 0;
        }
        //删除概率最大的节点
        int maxi = 0;
        double maxpro = -1.0;
        for (int i = 0; i < nnodes; i++)
        {
            if (pro[i] > maxpro)
            {
                maxpro = pro[i];
                maxi = i;
            }
        }
        outcore[maxi] = true;
        DEG[maxi] = 0;

        //他邻居的度数减少
        for (int e = head[maxi]; ~e; e = nxt[e])
            DEG[to[e]]--;
        clean();
    }
    for (int i = 0; i < nnodes; i++)
    {
        DEG[i] = tempdeg[i];
        outcore[i] = tempoutcore[i];
    }
    n = 0;
    for (int i = 0; i < nnodes; ++i) {
        if (!outcore[i])
            nV[i] = n, oriID[n++] = i;
    }
    //n为预处理后的节点数量
    //for(int i = 0; i < n; ++ i) Graph[i].reset();
    memset(head, -1, n * sizeof(int));
    memset(degree, 0, n * sizeof(int));
    memset(degree_temp, 0, n * sizeof(int));

    ecnt = 0;
    for (int i = 0; i < nedges; ++i) {
        int u = E[i].FI, v = E[i].SE;
        if (outcore[u] || outcore[v]) continue;
        u = nV[u]; v = nV[v];
        ++degree[u]; ++degree[v];
        ++degree_temp[u]; ++degree_temp[v];
        addedge(head + u, v);
        addedge(head + v, u);
        //Graph[u].set(v);
        //Graph[v].set(u);
    }

    /*
    for (int i = 0 ; i < n; i++)
        cout << degree[i] << endl;
    */
    cout << "预处理后的点数量为" << n << endl;
}

void nonpreprocess() {
    //Start to preprocess the graph    
    memset(DEG, 0, nnodes * sizeof(int));
    memset(head, -1, nnodes * sizeof(int));
    ecnt = 0;
    for (int i = 0; i < nedges; ++i) {
        ++DEG[E[i].FI];
        ++DEG[E[i].SE];
        addedge(head + E[i].FI, E[i].SE);
        addedge(head + E[i].SE, E[i].FI);
    }
    memset(outcore, 0, nnodes * sizeof(bool));
    n = 0;
    for (int i = 0; i < nnodes; ++i) {
        if (!outcore[i])
            nV[i] = n, oriID[n++] = i;
    }

    //for(int i = 0; i < n; ++ i) Graph[i].reset();
    memset(head, -1, n * sizeof(int));
    memset(degree, 0, n * sizeof(int));
    memset(degree_temp, 0, n * sizeof(int));

    ecnt = 0;
    for (int i = 0; i < nedges; ++i) {
        int u = E[i].FI, v = E[i].SE;
        if (outcore[u] || outcore[v]) continue;
        u = nV[u]; v = nV[v];
        ++degree[u]; ++degree[v];
        ++degree_temp[u]; ++degree_temp[v];
        addedge(head + u, v);
        addedge(head + v, u);
        //Graph[u].set(v);
        //Graph[v].set(u);
    }

}

void read(string filename) {
    ifstream inFile;
    inFile.open(filename);
    if (!inFile) {
        cout << "Unable to open file";
        exit(1); // terminate with error
    }
    //这里要求每个文件开头写上节点和边的数量
    inFile >> nnodes;
    inFile >> nedges;
    inFile >> startfromzero;
    //scanf("p edge %ld %ld", &nnodes, &nedges);
    cout << "网络节点数" << nnodes << " " << "网络边数" << nedges << endl;//
    for (int i = 0; i < nedges; ++i) {
        // int u, v; char e;
         //this is for the NDR
         //scanf("%c %d %d", &e, &E[i].FI, &E[i].SE);
        inFile >> E[i].FI;
        inFile >> E[i].SE;
        //this is for the DIMCAS
        //scanf("%d%d", &E[i].FI, &E[i].SE);
        if (!startfromzero)
        {
            --E[i].FI;
            --E[i].SE;
        }
    }
    inFile.close();
}

void read_weight(string filename)
{
    ifstream inFile;
    inFile.open(filename);
    if (!inFile) {
        cout << "Unable to open file";
        exit(1); // terminate with error
    }
    for (int i = 0; i < nnodes; i++)
        inFile >> vertex_weight[i];
    inFile.close();
}


void read_ratio(string filename, double ratio) {
    ifstream inFile;
    inFile.open(filename);
    if (!inFile) {
        cout << "Unable to open file";
        exit(1); // terminate with error
    }
    //这里要求每个文件开头写上节点和边的数量
    inFile >> nnodes;
    inFile >> nedges;
    inFile >> startfromzero;

    int ub_n = (int)nnodes * ratio, ub_e = 0;

    //scanf("p edge %ld %ld", &nnodes, &nedges);
    cout << "网络节点数" << nnodes << " " << "网络边数" << nedges << endl;//
    for (int i = 0; i < nedges; ++i) {

        inFile >> temp_E[i].FI;
        inFile >> temp_E[i].SE;
        if (temp_E[i].FI <= ub_n && temp_E[i].SE <= ub_n)
        {
            E[ub_e].FI = temp_E[i].FI;
            E[ub_e].SE = temp_E[i].SE;
            
            if (!startfromzero)
            {
                --E[ub_e].FI;
                --E[ub_e].SE;
            }
            ub_e++;
        }
    }
    nnodes = ub_n, nedges = ub_e;
    cout << ratio << "后的网络节点数" << nnodes << " " << "网络边数" << nedges << endl;//
    inFile.close();
}


double minw_heu(int q)
{
    double total_weight = 0;
    for (int i = 0; i < n; i++)
    {
        total_weight += vertex_weight[i];
    }
    printf("总权重为：%lf\n", total_weight);
    double upper_bound = 0;
    vector<int> temp_degree(n);
    for (int i = 0; i < n; i++)
        temp_degree[i] = degree[i];

    while (!outcore[q])
    {
        upper_bound = 0;
        for (int i = 0; i < n; i++)
        {
            if(!outcore[i])
            upper_bound += vertex_weight[i];
        }
        double max_weight = -1;
        int max_id = -1;
        for (int i = 0; i < n; i++)
        {
            if (vertex_weight[i] > max_weight && !outcore[i] && i != q)
            {
                max_weight = vertex_weight[i];
                max_id = i;
            }
        }

        if (max_id != -1)
        {
            outcore[max_id] = 1;
            for (int e = head[max_id]; ~e; e = nxt[e])
            {
                int v = to[e];
                if (!outcore[v]) temp_degree[v]--;
            }
            temp_degree[max_id] = 0;
            int fr = 0, re = 0;
            for (int i = 0; i < n; ++i) {
                if (temp_degree[i] < K) {
                    que[re++] = i;
                    outcore[i] = 1;//outcore[i]表示是否删除
                    /*exist_oldid[i] = 0;*/
                }
            }

            while (fr ^ re) {
                int u = que[fr++];
                for (int e = head[u]; ~e; e = nxt[e]) {
                    int v = to[e];
                    if (outcore[v]) continue;
                    if (--temp_degree[v] < K) {
                        que[re++] = v;
                        outcore[v] = true;
                        temp_degree[v] = 0;
                        /*exist_oldid[v] = 0;*/
                    }
                }
            }//不断重复删除
        }
    }
    memset(outcore, 0, n * sizeof(int));
    printf("权重的ub为：%lf\n", upper_bound);
    return upper_bound;
}

void addtoS_weight_enum(int u)
{
    ins[u] = 1;
    for (auto v : svex)
    {
        if (check_neighbor(u, v))
        {
            adj[u]++;
            adj[v]++;
        }
    }

}

void delfrS_weight_enum(int u)
{
    ins[u] = 0;
    adj[u] = 0;
    for (auto v : svex)
    {
        if (check_neighbor(u, v))
            adj[v]--;
    }
}

void delfrU_weight_enum(int u)
{
    del[u] = 1;
    degree[u] = 0;
    for (int e = head[u]; ~e; e = nxt[e])
    {
        int v = to[e];
        if (!del[v] && check_neighbor(u, v))
            degree[v]--;
    }
}

void addtoU_weight_enum(int u)
{
    del[u] = 0;
    for (int e = head[u]; ~e; e = nxt[e])
    {
        int v = to[e];
        if (!del[v] && check_neighbor(u, v))
        {
            degree[u]++;
            degree[v]++;
        }
    }
}

bool checkcore(vector<int> svex)
{
    for (auto u : svex)
    {
        if (!tempoutcore[u]) return true;
    }
    return false;
}

int minsize = MaxN;
void enum_weight(int q, int v, double& ub, clock_t start, bool add_flag, int &ans_size)
{
    if(v==-1) printf("s_size:%d\n", svex.size());
    if (add_flag)
    {
        cout << endl;
       /* printf("undel_list: ");
        for (int i = 0; i < n; i++)
        {
            if (!del[i]) cout << i << " ";
        }
        cout << endl;*/
      
        /*if (temp_weight >= ub)
        {
            cout << "overweight" << endl;
            return;
        }*/

        clock_t finish;
        finish = clock();
        if ((double)(finish - start) / CLOCKS_PER_SEC >= timelimit)
        {
            printf("Time exceeds\n");
            exit(100);
        }

        vector<int> temp_svex;
        temp_svex = svex;
        memset(tempadj, 0, sizeof(int) * n);
        memset(tempoutcore, 0, sizeof(int) * n);
        for (auto u : svex)
            tempadj[u] = adj[u];




        queue<int> que;
        for (auto u : temp_svex)
        {
            if (tempadj[u] < K)
                que.push(u);
        }

        while (!que.empty())
        {
            int u = que.front();
            que.pop();
            tempoutcore[u] = 1;
            for (auto w : temp_svex)
            {
                if (!tempoutcore[w] && check_neighbor(u, w))
                {
                    tempadj[w]--;
                    if (tempadj[w] < K) que.push(w);
                }
            }
        }

        if (checkcore(temp_svex) && !tempoutcore[q])
        {
            int cnt = 0;
            for (auto u : temp_svex)
            {
                if (!tempoutcore[u])
                {
                    cnt++;
                }
            }
            printf("size = %d\n", cnt);
            if (cnt < minsize) minsize = cnt;
            double temp_weight = 0;
            for (auto u : temp_svex)
            {
                if (!tempoutcore[u])
                    temp_weight += vertex_weight[u];
            }
            cout << "temp_weight = " << temp_weight << endl << endl;
            if (temp_weight < ub)
            {
                printf("Novel ub: %lf\n", temp_weight);
                ub = temp_weight;
                ans_size = cnt;
            }
        }
    }
 
    int u = v + 1;
    while (u < n)
    {
        if (u != q && !del[u]) break;
        u = u + 1;
    }
   
    if (u >= n) return;

  /*  if (v == q || del[v])
    {
        printf("%d has been deleted\n", v);
        enum_weight(q, v + 1, ub, start, 0, ans_size);
        return;
    }*/
    //分支1
   
    addtoS_weight_enum(u);
    svex.push_back(u);
    printf("add:current_v = %d,ub = %lf, init_size=%d, min_size = %d\n", u, ub, svex.size(), minsize);
    enum_weight(q, u, ub, start, 1, ans_size);
    svex.pop_back();
    delfrS_weight_enum(u);
    //删除U，并删除图中度数小于K的节点
    delfrU_weight_enum(u);
    printf("del:current_v = %d,ub = %lf, init_size=%d, min_size = %d\n", u, ub, svex.size(), minsize);


        //printf("s_size:%d\n", svex.size());
        //for (int i = 0; i < n; i++)
        //    cout << i << " " << degree[i] << endl;
    

    memset(degree_temp, 0, sizeof(int)* n);
    for (int i = 0; i < n; i++)
        degree_temp[i] = degree[i];

    vector<int> todel;
    queue<int> Q;
    for (int j = 0; j < n; j++)
    {
        if (!del[j] && degree_temp[j] < K)
        {
            if (ins[j])
            {
                printf("svex has been deleted, return\n");
                addtoU_weight_enum(u); //S中有节点度数小于K，剪枝
                return;
            }
            Q.push(j);
        }
    }

    while (!Q.empty())
    {
        int top = Q.front();
        Q.pop();
        todel.push_back(top);
        del[top] = 1;
        if (ins[top])
        {
           // if (u == 6) printf("degree 6 = %d\n", degree[6]);
            printf("svex is been deleted, return\n");
            for (auto x : todel) del[x] = 0;
            addtoU_weight_enum(u); //S中有节点被删除，剪枝
            //if (u == 6) printf("degree 6 = %d\n", degree[6]);
            return;
        }
        for (int j = 0; j < n; j++)
        {
            if (!del[j] && check_neighbor(top, j))
            {
                --degree_temp[j];
                if (degree_temp[j] < K)
                    Q.push(j);
            }
        }
    }
    for (int j = 0; j < todel.size(); j++)
        delfrU_weight_enum(todel[j]);
    //分支2
    enum_weight(q, u, ub, start, 0, ans_size);
    for (int j = todel.size() - 1; j >= 0; j--)
        addtoU_weight_enum(todel[j]);
    addtoU_weight_enum(u);
}

inline void delfrD(int u) {
    del[u] = true;
    delvex.push_back(u);
    for (int e = head[u]; ~e; e = nxt[e]) --degree[to[e]];
}

inline void addtoD(int u) {
    del[u] = false;
    delvex.pop_back();
    for (int e = head[u]; ~e; e = nxt[e]) ++degree[to[e]];
}

inline void delfrS(int u) {
    ins[u] = false;
    svex.pop_back();
    for (int i = 0; i < n; ++i) {
        if (i != u) {
            if (!checkedge(u, i))
                --notadj[i];
            else if (ins[i])
                --adj[i];
        }
    }
    adj[u] = 0;
}

inline void addtoS(int u) {
    ins[u] = true;
    svex.push_back(u);
    adj[u] = 0;
    for (int i = 0; i < n; ++i) {
        if (i != u) {
            if (!checkedge(u, i))
                ++notadj[i];
            else if (ins[i]) {
                ++adj[i];
                ++adj[u];
            }
        }
    }
}

inline bool canadd(int u) {
    int tot = 0;
    for (auto v : svex) if (!checkedge(u, v)) {
        if (++tot >= S) return false;
        if (notadj[v] >= S - 1) return false;
    }
    return true;
}


int dis[Maxn];

void bfs(int s) {
    memset(dis, -1, n * sizeof(int));
    dis[s] = 0;
    int fr = 0, re = 0;
    que[re++] = s;
    while (fr ^ re) {
        int u = que[fr++];
        for (int e = head[u]; ~e; e = nxt[e]) {
            int v = to[e];
            if (del[v] || dis[v] != -1) continue;
            dis[v] = dis[u] + 1;
            que[re++] = v;
        }
    }
}

bool bb_s_plex(int curS) {
    //cout << "The current S: " << S << endl;

    if (curS < LB) return false;
    int minID = -1;
    for (int i = 0; i < n; ++i) if (!del[i]) {
        if (minID == -1 || degree[i] < degree[minID])
            minID = i;
    }

    if (degree[minID] >= curS - S)
    {
        for (int i = 0; i < n; i++)
        {
            if (!del[i])
            {
                printf("其中之一的节点是：%d\n", oriID[i]);
                break;
            }
        }
        output_min_k_core();
        cout << "整个图都是splex，返回正确" << endl;
        cout << "此时图中有" << curS << "个节点" << endl;
        cout << "S中有" << svex.size() << "个节点" << endl;
        return true;
    }
    if (degree[minID] <= LB - 1 - S) {
        if (ins[minID])
            return false;
        delfrD(minID);
        bool ret = bb_s_plex(curS - 1);
        addtoD(minID);
        return ret;
    }

    int maxID = -1;
    for (int i = 0; i < n; ++i) if (!del[i]) {
        if (maxID == -1 || notadj[i] > notadj[maxID])
            maxID = i;
    }

    if (!ins[maxID] && notadj[maxID] >= S) {
        delfrD(maxID);
        bool ret = bb_s_plex(curS - 1);
        addtoD(maxID);
        return ret;
    }

    if (ins[maxID] && notadj[maxID] >= S - 1) {
        if (notadj[maxID] >= S) return false;
        std::vector <int> todel;
        for (int i = 0; i < n; ++i) if (!del[i] && !ins[i] && !checkedge(maxID, i)) {
            todel.push_back(i);
        }
        if (todel.size()) {
            for (auto x : todel) delfrD(x);
            bool ret = bb_s_plex(curS - todel.size());
            for (auto x : todel) addtoD(x);
            return ret;
        }
    }

    // This part (bfs) can be made recursively such that we can delete more useless vertices
    std::set <int> sofar;
    for (auto u : svex) {
        bfs(u);
        for (int v = 0; v < n; ++v) if (!del[v]) {
            if (dis[v] == -1 || dis[v] > fmax(2, S + S - LB + 1)) {
                if (ins[v]) return false;
                sofar.insert(v);
            }
        }
    }
    if (sofar.size()) {
        for (auto x : sofar) delfrD(x);
        bool ret = bb_s_plex(curS - sofar.size());
        for (auto x : sofar) addtoD(x);
        return ret;
    }
    std::vector <int> branch;
    for (int x = 0; x < n; ++x) if (!del[x] && x != minID && !ins[x]) {
        if (!checkedge(minID, x))
            branch.push_back(x);
    }

    random_shuffle(branch.begin(), branch.end());
    if (ins[minID]) {
        int canselect = S - 1 - notadj[minID], pos = -1;
        bool ret = false;
        for (int i = 0; !ret && i < canselect; ++i) {
            delfrD(branch[i]);
            if (i && !canadd(branch[i - 1])) {
                addtoD(branch[i]);
                break;
            }
            if (i) {
                addtoS(branch[i - 1]);
                pos = i - 1;
            }
            ret |= bb_s_plex(curS - 1);
            addtoD(branch[i]);
        }
        if (ret) {
            for (int i = 0; i <= pos; ++i) delfrS(branch[i]);
            return true;
        }
        for (int i = canselect; i < (int)branch.size(); ++i) {
            delfrD(branch[i]);
        }
        if (canselect == 0 || canadd(branch[canselect - 1])) {
            if (canselect) addtoS(branch[canselect - 1]);
            ret |= bb_s_plex(curS - branch.size() + canselect);
            if (canselect) delfrS(branch[canselect - 1]);
        }
        for (int i = 0; i <= pos; ++i) delfrS(branch[i]);
        for (int i = canselect; i < (int)branch.size(); ++i) {
            addtoD(branch[i]);
        }
        return ret;
    }
    else {
        delfrD(minID);

        bool ret = bb_s_plex(curS - 1);
        addtoD(minID);
        if (ret) return true;
        int canselect = S - 1 - notadj[minID];
        if (!canadd(minID)) return false;
        addtoS(minID);
        int pos = -1;
        for (int i = 0; !ret && i < canselect; ++i) {
            delfrD(branch[i]);
            if (i && !canadd(branch[i - 1])) {
                addtoD(branch[i]);
                break;
            }
            if (i) {
                addtoS(branch[i - 1]);
                pos = i - 1;
            }
            ret |= bb_s_plex(curS - 1);
            addtoD(branch[i]);
        }
        if (ret) {
            for (int i = 0; i <= pos; ++i) delfrS(branch[i]);
            delfrS(minID);
            return true;
        }
        for (int i = canselect; i < (int)branch.size(); ++i) {
            delfrD(branch[i]);
        }
        if (canselect == 0 || canadd(branch[canselect - 1])) {
            if (canselect) addtoS(branch[canselect - 1]);
            ret |= bb_s_plex(curS - branch.size() + canselect);
            if (canselect) delfrS(branch[canselect - 1]);
        }
        for (int i = 0; i <= pos; ++i) delfrS(branch[i]);
        for (int i = canselect; i < (int)branch.size(); ++i) {
            addtoD(branch[i]);
        }
        delfrS(minID);
        return ret;
    }
}


int iterative() {
    int y = n; S = 0;
    for (int i = 0; i < n; ++i) {
        S++;
        cout << "S: " << S << endl;
        LB = K + S;
        delvex.clear(); svex.clear();

        memset(ins, 0, n * sizeof(bool));
        memset(del, 0, n * sizeof(bool));
        memset(notadj, 0, n * sizeof(int));

        if (bb_s_plex(n)) return S + K;
        for (int i = 0; i < n; i++) {
            degree[i] = degree_temp[i];
        }


    }
    return 0;
}

int oldhead[MaxN], oldto[2 * MaxM], oldnxt[2 * MaxM], olddegree[MaxN], oldoutcore[MaxN], oldn;
double old_weight[MaxN];
int* temp_oldexist;
//memset(temp_oldexist, 1, sizeof(int)* nnodes);


void storeoldgraph()
{
    temp_oldexist = new int[nnodes];
    memset(temp_oldexist, 1, sizeof(int) * nnodes);
    oldn = n;
    for (int i = 0; i < nnodes; i++)
    {
        temp_oldexist[i] = exist_oldid[i];
        old_weight[i] = vertex_weight[i];
    }
    for (int i = 0; i < nnodes; i++)
    {
        oldhead[i] = head[i];
        olddegree[i] = degree[i];
        oldoutcore[i] = outcore[i];
    }
    for (int i = 0; i < 2 * nedges; i++)
    {
        oldto[i] = to[i];
        oldnxt[i] = nxt[i];
    }
}

void turn_to_oldgraph()
{
    n = oldn;
    for (int i = 0; i < nnodes; i++)
    {
        exist_oldid[i] = temp_oldexist[i];
        vertex_weight[i] = old_weight[i];
    }
    for (int i = 0; i < nnodes; i++)
    {
        head[i] = oldhead[i];
        degree[i] = olddegree[i];
        outcore[i] = oldoutcore[i];
    }
    for (int i = 0; i < 2 * nedges; i++)
    {
        to[i] = oldto[i];
        nxt[i] = oldnxt[i];
    }
}


void seek_remain_graph()
{
    int cnte = 0;
    for (int u = 0; u < n; u++)
    {
        if (!outcore[u])
        {
            for (int e = head[u]; ~e; e = nxt[e])
            {
                int v = to[e];
                if (!outcore[v] && v > u)
                {
                    newE[cnte].FI = u;
                    newE[cnte].SE = v;
                    cnte++;
                }
            }
        }
    }

    cout << "边数为" << cnte << endl;
    int tempn = 0;

    int* newid = new int[n];
    memset(newid, -1, n * sizeof(int));



    for (int i = 0; i < n; ++i) {
        if (!outcore[i])
        {
            newid[i] = tempn;
            vertex_weight[tempn] = vertex_weight[i];
            nV[oriID[i]] = tempn, oriID[tempn++] = oriID[i];// nV[i] = tempn
        }
        else exist_oldid[oriID[i]] = 0;
    }
    //n为预处理后的节点数量
    //for(int i = 0; i < n; ++ i) Graph[i].reset();
    memset(head, -1, tempn * sizeof(int));
    fill(nxt, nxt + 2 * MaxM, 0);
    fill(to, to + 2 * MaxM, 0);
    memset(degree, 0, tempn * sizeof(int));
    memset(degree_temp, 0, tempn * sizeof(int));

    ecnt = 0;
    for (int i = 0; i < cnte; ++i) {
        int u = newE[i].FI, v = newE[i].SE;
        if (outcore[u] || outcore[v]) continue;
        u = newid[u]; v = newid[v];
        ++degree[u]; ++degree[v];
        ++degree_temp[u]; ++degree_temp[v];

        addedge(head + u, v);
        addedge(head + v, u);
        //Graph[u].set(v);
        //Graph[v].set(u);
    }

    /*
    for (int j = 0; j < tempn; j++)
        cout << degree[j] << endl;
    */


    memset(outcore, 0, tempn * sizeof(bool));
    n = tempn;
    cout << "网络一开始的节点数量为" << nnodes << endl;
    for (int i = n; i < nnodes; i++)
        outcore[i] = 1;
    for (int i = 0; i < n; i++) nV[oriID[i]] = i;
    cout << "2nd_order处理后的点数量为" << n << endl;
}

int visit[MaxN];

int countintersect(int u, int v)
{
    memset(visit, 0, sizeof(int) * n);
    vector<int> ans;
    for (int e = head[u]; ~e; e = nxt[e])
    {
        if (!visit[to[e]])
        {
            visit[to[e]] = 1;
        }
    }
    for (int e = head[v]; ~e; e = nxt[e])
    {
        if (visit[to[e]])
        {
            ans.push_back(to[e]);
        }
    }
    /*set<int> uset, vset, ans;
    for (int e = head[u]; ~e; e = nxt[e])
        uset.insert(to[e]);
    for (int e = head[v]; ~e; e = nxt[e])
        vset.insert(to[e]);
    for (auto it = uset.begin(); it != uset.end(); it++)
    {
        if (vset.find(*it) != vset.end())
            ans.insert(*it);
    }*/
    return ans.size();
}
void countingtriangles(int* triangles)
{
    for (int u = 0; u < n; u++)
    {
        if (!outcore[u])
        {
            for (int e = head[u]; ~e; e = nxt[e])
            {
                int v = to[e];
                if (u < v)
                {
                    triangles[e] = countintersect(u, v);
                }
            }
        }
    }
}

void removevertex(int u)
{
    //cout << "remove" << endl;
    vector<int> neighbour;
    int laste;
    for (int e = head[u]; ~e; e = nxt[e])
        neighbour.push_back(to[e]);
    for (int i = 0; i < neighbour.size(); i++)
    {
        int v = neighbour[i];
        for (int e = head[v]; ~e; e = nxt[e])
        {
            if (to[e] == u)
            {
                if (e == head[v]) head[v] = nxt[e];
                else
                {
                    nxt[laste] = nxt[e];
                    to[e] = -1;
                }
                break;
            }
            laste = e;
        }
    }
    // degree[u] = 0;
    outcore[u] = true;
    head[u] = -1;
}

void removedges(int u, int v)
{
    int laste;
    for (int e = head[u]; ~e; e = nxt[e])
    {
        if (to[e] == v)
        {
            if (e == head[u]) head[u] = nxt[e];
            else
            {
                nxt[laste] = nxt[e];
                to[e] = -1;
            }
            break;
        }
        laste = e;
    }

    for (int e = head[v]; ~e; e = nxt[e])
    {
        if (to[e] == u)
        {
            if (e == head[v]) head[v] = nxt[e];
            else
            {
                nxt[laste] = nxt[e];
                to[e] = -1;
            }
            break;
        }
        laste = e;
    }
    degree[u]--;
    degree[v]--;
}

void second_order_reduction(int s, int lb)
{
    queue<int> q;
    for (int i = 0; i < n; i++)
    {
        if (degree[i] < lb - s && !outcore[i])
        {
            //cout << "entry queue!" << endl;
            q.push(i);
        }
    }
    while (true)
    {
        if (!q.empty())
        {
            //cout << "not empty" << endl;
            int u = q.front();
            q.pop();
            for (int e = head[u]; ~e; e = nxt[e])
            {
                //先减少邻居的度数
                int v = to[e];
                degree[v]--;
                if (degree[v] < lb - s && !outcore[v])
                    q.push(v);
            }
            removevertex(u);
        }
        else
        {
            fill(triangles, triangles + 2 * MaxM, 0);
            countingtriangles(triangles);

            set<pair<int, int>> edges;
            for (int u = 0; u < n; u++)
            {
                if (!outcore[u])
                {
                    for (int e = head[u]; ~e; e = nxt[e])
                    {
                        int v = to[e];
                        //cout << "triangle e is" << triangles[e] << endl;
                        if (u < v && triangles[e] + 2 * s < lb)
                        {
                            //cout << "remove edge" << endl;
                            //cout << u << " " << v << endl;
                            edges.insert(make_pair(u, v));
                        }
                    }
                }
            }

        /*    vector<int> vts;
            appear = new int[n];
            memset(appear, 0, sizeof(int) * n);*/
            
            set<int> vts;
            for (auto e : edges)
            {
                vts.insert(e.first);
                vts.insert(e.second);
                //if (!appear[e.first])
                //{
                //    vts.push_back(e.first);
                //    appear[e.first] = 1;
                //}

                //if (!appear[e.second])
                //{
                //    vts.push_back(e.second);
                //    appear[e.first] = 1;
                //}

                removedges(e.first, e.second);
            }
            for (auto v : vts)
            {
                if (degree[v] + s < lb)
                    q.push(v);
            }
            if (q.empty())
                break;
        }
    }
}

bool bb_s_plex_new(int curS, clock_t start, clock_t finish, bool& timeflag) {
    //cout << "The current S: " << S << endl;
    //s:解集 c:候选集 S:S-plex的S
    //notadj[i]：i在s中非邻居的数量
    //curS指的不是S大小，而是s+c的大小
    //LB：图中最少的节点个数（S+K)

    //没有用到论文2里面的预处理


   finish = clock();
   int time = (double)(finish - start) / CLOCKS_PER_SEC;
    
    if (time > 1800)
    {
        timeflag = false;
        return false;

    }
    

    //这里就可以返回了，没必要往后面继续跑了

    if (svex.size() >= LB)
    {
        //加上这三句试一试
        //finish = clock();
        finaltime = (double)(finish - start) / CLOCKS_PER_SEC;
        cout << "规模大于LB，返回正确，时间是：" << finaltime << endl;
        return true;
    }


    if (curS < LB)
    {
        return false;
    }
    int minID = -1;
    for (int i = 0; i < n; ++i) if (!del[i]) {
        if (minID == -1 || degree[i] < degree[minID])
            minID = i;
    }

    if (degree[minID] >= curS - S)
    {
        for (int i = 0; i < n; i++)
        {
            if (!del[i])
            {
                printf("其中之一的节点是：%d\n", oriID[i]);
                break;
            }
        }
        output_min_k_core();
        cout << "整个图都是splex，返回正确" << endl;
        cout << "此时图中有" << curS << "个节点" << endl;
        cout << "S中有" << svex.size() << "个节点" << endl;
        //output_min_k_core();
        return true; //说明整个图都是s-plex
    }
    if (degree[minID] <= LB - 1 - S) //说明MinID不可能加入最终解，第9行
    {
        if (ins[minID])
            return false;
        delfrD(minID);
        bool ret = bb_s_plex_new(curS - 1, start, finish, timeflag);
        addtoD(minID);
        return ret;
    }

    int maxID = -1;//maxID指在s中非邻居数量最多的节点
    for (int i = 0; i < n; ++i) if (!del[i]) {
        if (maxID == -1 || notadj[i] > notadj[maxID])
            maxID = i;
    }

    //如果maxID不在s中，且他在s中非邻居大于S，那他必被删，第5行，lemma2

    if (!ins[maxID] && notadj[maxID] >= S) {
        delfrD(maxID);
        bool ret = bb_s_plex_new(curS - 1, start, finish, timeflag);
        addtoD(maxID);
        return ret;
    }

    //如果maxID在s中，且非邻居大于S-1，该节点是临界点，可以借助临界点来删除c中和临界点不相邻的节点，第5行，lemma2

    if (ins[maxID] && notadj[maxID] >= S - 1) {
        if (notadj[maxID] >= S) return false;
        std::vector <int> todel;
        for (int i = 0; i < n; ++i) if (!del[i] && !ins[i] && !checkedge(maxID, i)) {
            todel.push_back(i);
        }
        if (todel.size()) {
            for (auto x : todel) delfrD(x);
            bool ret = bb_s_plex_new(curS - todel.size(), start, finish, timeflag);
            for (auto x : todel) addtoD(x);
            return ret;
        }
    }

    // This part (bfs) can be made recursively such that we can delete more useless vertices
    //c中和s中节点过远的节点也会被删除，第11行
    std::set <int> sofar;
    for (auto u : svex) {
        bfs(u);
        for (int v = 0; v < n; ++v) if (!del[v]) {
            if (dis[v] == -1 || dis[v] > fmax(2, S + S - LB + 1)) {
                if (ins[v]) return false;
                sofar.insert(v);
            }
        }
    }
    if (sofar.size()) {
        for (auto x : sofar) delfrD(x);
        bool ret = bb_s_plex_new(curS - sofar.size(), start, finish, timeflag);
        for (auto x : sofar) addtoD(x);
        return ret;
    }

    //开始分支（以minID为基准，在这里指的是k-unsatisfied节点，因为若整个图不满足780行的情况，那么minID肯定是k-unsatisfied)
    std::vector <int> branch;
    for (int x = 0; x < n; ++x) if (!del[x] && x != minID && !ins[x]) {
        if (!checkedge(minID, x))
            branch.push_back(x);
    }

    random_shuffle(branch.begin(), branch.end());
    //若minID在s中
    if (ins[minID]) {
        int canselect = S - 1 - notadj[minID], pos = -1;
        bool ret = false;
        for (int i = 0; !ret && i < canselect; ++i) {
            delfrD(branch[i]);
            if (i && !canadd(branch[i - 1])) {
                addtoD(branch[i]);
                break;
            }
            if (i) {
                addtoS(branch[i - 1]);
                pos = i - 1;
            }
            ret |= bb_s_plex_new(curS - 1, start, finish, timeflag);
            addtoD(branch[i]);
        }
        if (ret) {
            for (int i = 0; i <= pos; ++i) delfrS(branch[i]);
            return true;
        }
        //后续与minID不相连的节点都得删除
        for (int i = canselect; i < (int)branch.size(); ++i) {
            delfrD(branch[i]);
        }
        //最后一个分支
        if (canselect == 0 || canadd(branch[canselect - 1])) {
            if (canselect) addtoS(branch[canselect - 1]);
            ret |= bb_s_plex_new(curS - branch.size() + canselect, start, finish, timeflag);
            if (canselect) delfrS(branch[canselect - 1]);
        }
        for (int i = 0; i <= pos; ++i) delfrS(branch[i]);
        for (int i = canselect; i < (int)branch.size(); ++i) {
            addtoD(branch[i]);
        }
        return ret;
    }

    //若minID不在s中
    else {
        delfrD(minID);

        bool ret = bb_s_plex_new(curS - 1, start, finish, timeflag);
        addtoD(minID);

        if (ret) return true;
        int canselect = S - 1 - notadj[minID];
        if (!canadd(minID)) return false;
        addtoS(minID);
        int pos = -1;
        for (int i = 0; !ret && i < canselect; ++i) {
            delfrD(branch[i]);
            if (i && !canadd(branch[i - 1])) {
                addtoD(branch[i]);
                break;
            }
            if (i) {
                addtoS(branch[i - 1]);
                pos = i - 1;
            }
            ret |= bb_s_plex_new(curS - 1, start, finish, timeflag);
            addtoD(branch[i]);
        }
        if (ret) {
            for (int i = 0; i <= pos; ++i) delfrS(branch[i]);
            delfrS(minID);
            return true;
        }
        for (int i = canselect; i < (int)branch.size(); ++i) {
            delfrD(branch[i]);
        }
        if (canselect == 0 || canadd(branch[canselect - 1])) {
            if (canselect) addtoS(branch[canselect - 1]);
            ret |= bb_s_plex_new(curS - branch.size() + canselect, start, finish, timeflag);
            if (canselect) delfrS(branch[canselect - 1]);
        }
        for (int i = 0; i <= pos; ++i) delfrS(branch[i]);
        for (int i = canselect; i < (int)branch.size(); ++i) {
            addtoD(branch[i]);
        }
        delfrS(minID);
        return ret;
    }
}

long long factorial(long long n)
{
    long long result = 1;
    for (long long i = 2; i <= n; i++) {
        result *= i;
    }
    return result;
}

long long combination(long long n, long long r)
{
    if (r > n) return 0;
    if (r == n || r == 0) return 1;
    return factorial(n) / (factorial(r) * factorial(n - r));
}

vector<vector<int>> topr_ans;

void backtrack(vector<int>& nums, vector<int>& cur, int r, int start, int target_num)
{
    if (topr_ans.size() >= target_num) return;

    if (cur.size() == r)
    {
        for (int i = 0; i < topr_ans.size(); i++)
        {
            if (cur == topr_ans[i])
                return;
        }
        topr_ans.push_back(cur);
        
        return;
    }

    for (int i = start; i < nums.size(); i++)
    {
        cur.push_back(nums[i]);
        backtrack(nums, cur, r, i + 1,target_num);
        cur.pop_back();
    }
}

bool cmp_weight(int a,int b)
{
    return vertex_weight[a] < vertex_weight[b];
}



bool bb_s_plex_topr(int curS, clock_t start, clock_t finish, bool& timeflag, int r, long long & c) {
    //cout << "The current S: " << S << endl;
    //s:解集 c:候选集 S:S-plex的S
    //notadj[i]：i在s中非邻居的数量
    //curS指的不是S大小，而是s+c的大小
    //LB：图中最少的节点个数（S+K)

    //没有用到论文2里面的预处理


  // finish = clock();
   //int time = (double)(finish - start) / CLOCKS_PER_SEC;
    /*
    if (time > 3600)
    {
        timeflag = false;
       return false;

    }
    */

    //这里就可以返回了，没必要往后面继续跑了

    if (svex.size() >= LB)
    {
        //加上这三句试一试
        finish = clock();
        finaltime = (double)(finish - start) / CLOCKS_PER_SEC;
        cout << "规模大于LB，返回正确，时间是：" << finaltime << endl;

        vector<int> nums, cur;
        for (int i = 0; i < svex.size(); i++)
        {
            if (!del[svex[i]])
                nums.push_back(svex[i]);
        }
        backtrack(nums, cur, LB, 0, r);
        c = topr_ans.size();
        printf("Find %d ans! Size is %d! Remain %d\n\n", c, LB, (r - c) > 0 ? r - c : 0);
        printf("用时: %lf\n\n", finaltime);
        if (r == 1)
        {
            vector<int> ans_temp;
            ans_temp = topr_ans[0];
            sort(ans_temp.begin(), ans_temp.end());
            ofstream outfile;
            outfile.open("/home/ljw/zqf/dataset/min_k_core.txt");
            for (int ii = 0; ii < ans_temp.size(); ii++)
            {
                for (int jj = ii + 1; jj < ans_temp.size(); jj++)
                {
                    if (check_neighbor(ans_temp[ii], ans_temp[jj]))
                        outfile << oriID[ans_temp[ii]] << " " << oriID[ans_temp[jj]] << endl;
                }
            }
        }
        if (c >= r)return true;
        return false;
    }


    if (curS < LB)
    {
        return false;
    }
    int minID = -1;
    for (int i = 0; i < n; ++i) if (!del[i]) {
        if (minID == -1 || degree[i] < degree[minID])
            minID = i;
    }

    if (degree[minID] >= curS - S)
    {
        cout << "整个图都是splex，返回正确" << endl;
        cout << "此时图中有" << curS << "个节点" << endl;
        cout << "S中有" << svex.size() << "个节点" << endl;
        //output_min_k_core();
        //int x = curS - 1, y = LB - 1;
        int x = curS, y = LB;  //不是query，没必要-1

        vector<int> nums, cur;
        for (int i = 0; i < n; i++)
        {
            if (!del[i])
                nums.push_back(i);
        }
       /* combination(x, y);*/
 /*       if (combination(x, y) <= 0 || combination(x, y) >= r)
        {
            printf("n=%d,r=%d, comb=%d,c=%d\n\n", x, y, combination(x, y), c);
            printf("find ans!\n");
            backtrack(nums, cur, y, 0, 1, r);
            return true;
        }*/
       // c += combination(x, y);
        backtrack(nums, cur, y, 0, r);
        c = topr_ans.size();
        finish = clock();
        double temp_time = (double)(finish - start) / CLOCKS_PER_SEC;
        printf("Already find %d ans! Size is %d! Remain %d\n", c, LB, (r - c) > 0 ? r - c : 0);
        printf("用时: %lf\n\n", temp_time);

        if (r == 1)
        {
            vector<int> ans_temp;
            ans_temp = topr_ans[0];
            sort(ans_temp.begin(), ans_temp.end());
            ofstream outfile;
            outfile.open("/home/ljw/zqf/dataset/min_k_core.txt");
            for (int ii = 0; ii < ans_temp.size(); ii++)
            {
                for (int jj = ii + 1; jj < ans_temp.size(); jj++)
                {
                    if (check_neighbor(ans_temp[ii], ans_temp[jj]) && startfromzero)
                        outfile << oriID[ans_temp[ii]] << " " << oriID[ans_temp[jj]] << endl;
                    else if (check_neighbor(ans_temp[ii], ans_temp[jj]) && !startfromzero)
                        outfile << oriID[ans_temp[ii]] + 1 << " " << oriID[ans_temp[jj]] + 1 << endl;
                }
            }
            /*for (int ii = 0; ii < ans_temp.size(); ii++)
            {
                for (int jj = ii + 1; jj < ans_temp.size(); jj++)
                {
                    if (check_neighbor(ans_temp[ii], ans_temp[jj]))
                        outfile << ans_temp[ii] << " " << ans_temp[jj] << endl;
                }
            }*/
        }

        if (c >= r) return true;
        return false;
       // return true; //说明整个图都是s-plex
    }
    if (degree[minID] <= LB - 1 - S) //说明MinID不可能加入最终解，第9行
    {
        if (ins[minID])
            return false;
        delfrD(minID);
        bool ret = bb_s_plex_topr(curS - 1, start, finish, timeflag, r, c);
        addtoD(minID);
        return ret;
    }

    int maxID = -1;//maxID指在s中非邻居数量最多的节点
    for (int i = 0; i < n; ++i) if (!del[i]) {
        if (maxID == -1 || notadj[i] > notadj[maxID])
            maxID = i;
    }

    //如果maxID不在s中，且他在s中非邻居大于S，那他必被删，第5行，lemma2

    if (!ins[maxID] && notadj[maxID] >= S) {
        delfrD(maxID);
        bool ret = bb_s_plex_topr(curS - 1, start, finish, timeflag, r, c);
        addtoD(maxID);
        return ret;
    }

    //如果maxID在s中，且非邻居大于S-1，该节点是临界点，可以借助临界点来删除c中和临界点不相邻的节点，第5行，lemma2

    if (ins[maxID] && notadj[maxID] >= S - 1) {
        if (notadj[maxID] >= S) return false;
        std::vector <int> todel;
        for (int i = 0; i < n; ++i) if (!del[i] && !ins[i] && !checkedge(maxID, i)) {
            todel.push_back(i);
        }
        if (todel.size()) {
            for (auto x : todel) delfrD(x);
            bool ret = bb_s_plex_topr(curS - todel.size(), start, finish, timeflag, r, c);
            for (auto x : todel) addtoD(x);
            return ret;
        }
    }

    // This part (bfs) can be made recursively such that we can delete more useless vertices
    //c中和s中节点过远的节点也会被删除，第11行
    std::set <int> sofar;
    for (auto u : svex) {
        bfs(u);
        for (int v = 0; v < n; ++v) if (!del[v]) {
            if (dis[v] == -1 || dis[v] > fmax(2, S + S - LB + 1)) {
                if (ins[v]) return false;
                sofar.insert(v);
            }
        }
    }
    if (sofar.size()) {
        for (auto x : sofar) delfrD(x);
        bool ret = bb_s_plex_topr(curS - sofar.size(), start, finish, timeflag, r, c);
        for (auto x : sofar) addtoD(x);
        return ret;
    }

    //开始分支（以minID为基准，在这里指的是k-unsatisfied节点，因为若整个图不满足780行的情况，那么minID肯定是k-unsatisfied)
    std::vector <int> branch;
    for (int x = 0; x < n; ++x) if (!del[x] && x != minID && !ins[x]) {
        if (!checkedge(minID, x))
            branch.push_back(x);
    }

    random_shuffle(branch.begin(), branch.end());
    //若minID在s中
    if (ins[minID]) {
        int canselect = S - 1 - notadj[minID], pos = -1;
        bool ret = false;
        for (int i = 0; !ret && i < canselect; ++i) {
            delfrD(branch[i]);
            if (i && !canadd(branch[i - 1])) {
                addtoD(branch[i]);
                break;
            }
            if (i) {
                addtoS(branch[i - 1]);
                pos = i - 1;
            }
            ret |= bb_s_plex_topr(curS - 1, start, finish, timeflag, r, c);
            addtoD(branch[i]);
        }
        if (ret) {
            for (int i = 0; i <= pos; ++i) delfrS(branch[i]);
            return true;
        }
        //后续与minID不相连的节点都得删除
        for (int i = canselect; i < (int)branch.size(); ++i) {
            delfrD(branch[i]);
        }
        //最后一个分支
        if (canselect == 0 || canadd(branch[canselect - 1])) {
            if (canselect) addtoS(branch[canselect - 1]);
            ret |= bb_s_plex_topr(curS - branch.size() + canselect, start, finish, timeflag, r, c);
            if (canselect) delfrS(branch[canselect - 1]);
        }
        for (int i = 0; i <= pos; ++i) delfrS(branch[i]);
        for (int i = canselect; i < (int)branch.size(); ++i) {
            addtoD(branch[i]);
        }
        return ret;
    }

    //若minID不在s中
    else {
        delfrD(minID);

        bool ret = bb_s_plex_topr(curS - 1, start, finish, timeflag, r, c);
        addtoD(minID);

        if (ret) return true;
        int canselect = S - 1 - notadj[minID];
        if (!canadd(minID)) return false;
        addtoS(minID);
        int pos = -1;
        for (int i = 0; !ret && i < canselect; ++i) {
            delfrD(branch[i]);
            if (i && !canadd(branch[i - 1])) {
                addtoD(branch[i]);
                break;
            }
            if (i) {
                addtoS(branch[i - 1]);
                pos = i - 1;
            }
            ret |= bb_s_plex_topr(curS - 1, start, finish, timeflag, r, c);
            addtoD(branch[i]);
        }
        if (ret) {
            for (int i = 0; i <= pos; ++i) delfrS(branch[i]);
            delfrS(minID);
            return true;
        }
        for (int i = canselect; i < (int)branch.size(); ++i) {
            delfrD(branch[i]);
        }
        if (canselect == 0 || canadd(branch[canselect - 1])) {
            if (canselect) addtoS(branch[canselect - 1]);
            ret |= bb_s_plex_topr(curS - branch.size() + canselect, start, finish, timeflag, r, c);
            if (canselect) delfrS(branch[canselect - 1]);
        }
        for (int i = 0; i <= pos; ++i) delfrS(branch[i]);
        for (int i = canselect; i < (int)branch.size(); ++i) {
            addtoD(branch[i]);
        }
        delfrS(minID);
        return ret;
    }
}

bool bb_s_plex_weight(int curS, clock_t start, clock_t finish, bool& timeflag, bool& flag, double &ub, int q) {


   finish = clock();
   int time = (double)(finish - start) / CLOCKS_PER_SEC;
   if (time > timelimit)
   {
       timeflag = false;
       return false;
   }


    if (svex.size() >= LB)
    {
        //加上这三句试一试
        finish = clock();
        finaltime = (double)(finish - start) / CLOCKS_PER_SEC;
        cout << "规模大于LB，返回正确，时间是：" << finaltime << endl;

        flag = true;
        sort(svex.begin(), svex.end(), cmp_weight);

        bool isq = false;
        double total_weight = 0;
        for (int i = 0; i < LB - 1; i++)
        {
            if (svex[i] == q) isq = true;
            total_weight += vertex_weight[svex[i]];
        }

        if (isq) total_weight += vertex_weight[svex[LB - 1]];
        else total_weight += vertex_weight[q];
        if (total_weight < ub) ub = total_weight;       
        return false;
    }


    if (curS < LB)
    {
        return false;
    }
    int minID = -1;
    for (int i = 0; i < n; ++i) if (!del[i]) {
        if (minID == -1 || degree[i] < degree[minID])
            minID = i;
    }

    if (degree[minID] >= curS - S)
    {
        cout << "整个图都是splex，返回正确" << endl;
        cout << "此时图中有" << curS << "个节点" << endl;
        cout << "S中有" << svex.size() << "个节点" << endl;

        flag = true;
        vector<int> ans_vex;
        for (int i = 0; i < n; i++)
        {
            if (!del[i])
                ans_vex.push_back(i);
        }
        sort(ans_vex.begin(), ans_vex.end(), cmp_weight);

        bool isq = false;
        double total_weight = 0;
        for (int i = 0; i < LB - 1; i++)
        {
            if (ans_vex[i] == q) isq = true;
            total_weight += vertex_weight[ans_vex[i]];
        }
        if (isq) total_weight += vertex_weight[ans_vex[LB - 1]];
        else total_weight += vertex_weight[q];
        if (total_weight < ub) ub = total_weight;

        return false;
        // return true; //说明整个图都是s-plex
    }
    if (degree[minID] <= LB - 1 - S) //说明MinID不可能加入最终解，第9行
    {
        if (ins[minID])
            return false;
        delfrD(minID);
        bool ret = bb_s_plex_weight(curS - 1, start, finish, timeflag, flag, ub, q);
        addtoD(minID);
        return ret;
    }

    int maxID = -1;//maxID指在s中非邻居数量最多的节点
    for (int i = 0; i < n; ++i) if (!del[i]) {
        if (maxID == -1 || notadj[i] > notadj[maxID])
            maxID = i;
    }

    //如果maxID不在s中，且他在s中非邻居大于S，那他必被删，第5行，lemma2

    if (!ins[maxID] && notadj[maxID] >= S) {
        delfrD(maxID);
        bool ret = bb_s_plex_weight(curS - 1, start, finish, timeflag, flag, ub, q);
        addtoD(maxID);
        return ret;
    }

    //如果maxID在s中，且非邻居大于S-1，该节点是临界点，可以借助临界点来删除c中和临界点不相邻的节点，第5行，lemma2

    if (ins[maxID] && notadj[maxID] >= S - 1) {
        if (notadj[maxID] >= S) return false;
        std::vector <int> todel;
        for (int i = 0; i < n; ++i) if (!del[i] && !ins[i] && !checkedge(maxID, i)) {
            todel.push_back(i);
        }
        if (todel.size()) {
            for (auto x : todel) delfrD(x);
            bool ret = bb_s_plex_weight(curS - todel.size(), start, finish, timeflag, flag, ub, q);
            for (auto x : todel) addtoD(x);
            return ret;
        }
    }

    // This part (bfs) can be made recursively such that we can delete more useless vertices
    //c中和s中节点过远的节点也会被删除，第11行
    std::set <int> sofar;
    for (auto u : svex) {
        bfs(u);
        for (int v = 0; v < n; ++v) if (!del[v]) {
            if (dis[v] == -1 || dis[v] > fmax(2, S + S - LB + 1)) {
                if (ins[v]) return false;
                sofar.insert(v);
            }
        }
    }
    if (sofar.size()) {
        for (auto x : sofar) delfrD(x);
        bool ret = bb_s_plex_weight(curS - sofar.size(), start, finish, timeflag, flag, ub, q);
        for (auto x : sofar) addtoD(x);
        return ret;
    }

    //开始分支（以minID为基准，在这里指的是k-unsatisfied节点，因为若整个图不满足780行的情况，那么minID肯定是k-unsatisfied)
    std::vector <int> branch;
    for (int x = 0; x < n; ++x) if (!del[x] && x != minID && !ins[x]) {
        if (!checkedge(minID, x))
            branch.push_back(x);
    }

    random_shuffle(branch.begin(), branch.end());
    //若minID在s中
    if (ins[minID]) {
        int canselect = S - 1 - notadj[minID], pos = -1;
        bool ret = false;
        for (int i = 0; !ret && i < canselect; ++i) {
            delfrD(branch[i]);
            if (i && !canadd(branch[i - 1])) {
                addtoD(branch[i]);
                break;
            }
            if (i) {
                addtoS(branch[i - 1]);
                pos = i - 1;
            }
            ret |= bb_s_plex_weight(curS - 1, start, finish, timeflag, flag, ub, q);
            addtoD(branch[i]);
        }
        if (ret) {
            for (int i = 0; i <= pos; ++i) delfrS(branch[i]);
            return true;
        }
        //后续与minID不相连的节点都得删除
        for (int i = canselect; i < (int)branch.size(); ++i) {
            delfrD(branch[i]);
        }
        //最后一个分支
        if (canselect == 0 || canadd(branch[canselect - 1])) {
            if (canselect) addtoS(branch[canselect - 1]);
            ret |= bb_s_plex_weight(curS - branch.size() + canselect, start, finish, timeflag, flag, ub, q);
            if (canselect) delfrS(branch[canselect - 1]);
        }
        for (int i = 0; i <= pos; ++i) delfrS(branch[i]);
        for (int i = canselect; i < (int)branch.size(); ++i) {
            addtoD(branch[i]);
        }
        return ret;
    }

    //若minID不在s中
    else {
        delfrD(minID);

        bool ret = bb_s_plex_weight(curS - 1, start, finish, timeflag, flag, ub, q);
        addtoD(minID);

        if (ret) return true;
        int canselect = S - 1 - notadj[minID];
        if (!canadd(minID)) return false;
        addtoS(minID);
        int pos = -1;
        for (int i = 0; !ret && i < canselect; ++i) {
            delfrD(branch[i]);
            if (i && !canadd(branch[i - 1])) {
                addtoD(branch[i]);
                break;
            }
            if (i) {
                addtoS(branch[i - 1]);
                pos = i - 1;
            }
            ret |= bb_s_plex_weight(curS - 1, start, finish, timeflag, flag, ub, q);
            addtoD(branch[i]);
        }
        if (ret) {
            for (int i = 0; i <= pos; ++i) delfrS(branch[i]);
            delfrS(minID);
            return true;
        }
        for (int i = canselect; i < (int)branch.size(); ++i) {
            delfrD(branch[i]);
        }
        if (canselect == 0 || canadd(branch[canselect - 1])) {
            if (canselect) addtoS(branch[canselect - 1]);
            ret |= bb_s_plex_weight(curS - branch.size() + canselect, start, finish, timeflag, flag, ub, q);
            if (canselect) delfrS(branch[canselect - 1]);
        }
        for (int i = 0; i <= pos; ++i) delfrS(branch[i]);
        for (int i = canselect; i < (int)branch.size(); ++i) {
            addtoD(branch[i]);
        }
        delfrS(minID);
        return ret;
    }
}


int best_iterative() {
    clock_t start, finish = 0;
    start = clock();
    bool timeflag = true;
    int y = n; S = 0;
    storeoldgraph();

    bool* tempnxt = new bool[n];
    memset(tempnxt, 0, n * sizeof(bool));
    int temp_n = n;
    //这里是不是应该改成n-k?
    for (int i = 0; i < temp_n - K; ++i) {
        S++;
        cout << "S: " << S << endl;
        cout << "S + K: " << S + K << endl;
        LB = K + S;

        memset(outcore, 0, n * sizeof(bool));
        second_order_reduction(S, LB);
        seek_remain_graph();



        delvex.clear(); svex.clear();


        memset(ins, 0, n * sizeof(bool));
        memset(del, 0, n * sizeof(bool));
        memset(notadj, 0, n * sizeof(int));
        if (bb_s_plex_new(n, start, finish, timeflag))
        {
            if (timeflag == false)
                return 0;
            break;
        }

        for (int i = 0; i < temp_n; i++) {
            degree[i] = degree_temp[i];
        }
        turn_to_oldgraph();
    }
    return K + S;
}

int best_iterative_all() {
    clock_t start, finish = 0;
    start = clock();
    bool timeflag = true;
    int y = n; S = 0;
    int ans = nnodes;
    storeoldgraph();

    bool* tempnxt = new bool[n];
    memset(tempnxt, 0, n * sizeof(bool));
    //这里是不是应该改成n-k?
    int temp_n = n;
    for (int i = 0; i < temp_n - K; ++i) {
        S++;
        cout << "S: " << S << endl;
        cout << "S + K: " << S + K << endl;
        LB = K + S;

        memset(outcore, 0, n * sizeof(bool));
        second_order_reduction(S, LB);
        seek_remain_graph();

        finish = clock();
        int time = (double)(finish - start) / CLOCKS_PER_SEC;

        if (time > 10000)
        {
            timeflag = false;
            return -1;

        }

        delvex.clear(); svex.clear();


        memset(ins, 0, n * sizeof(bool));
        memset(del, 0, n * sizeof(bool));
        memset(notadj, 0, n * sizeof(int));
        if (bb_s_plex_new(n, start, finish, timeflag))
        {
            if (LB < ans) ans = LB;
            if (timeflag == false)
                return 0;
        }

        for (int i = 0; i < temp_n; i++) {
            degree[i] = degree_temp[i];
        }
        turn_to_oldgraph();
    }
    return ans;
}

int best_iterative_all_query(int q) {
    clock_t start, finish = 0;
    start = clock();
    bool timeflag = true;
    int y = n; S = 0;
    int ans = nnodes;
    storeoldgraph();

    bool* tempnxt = new bool[n];
    memset(tempnxt, 0, n * sizeof(bool));
    int temp_n = n;
    for (int i = 0; i < temp_n - K; ++i) {
        S++;
        cout << "S: " << S << endl;
        cout << "S + K: " << S + K << endl;
        LB = K + S;

        memset(outcore, 0, n * sizeof(bool));
        second_order_reduction(S, LB);
        seek_remain_graph();



        delvex.clear(); svex.clear();
        if (!exist_oldid[q])
        {
            printf("Query vertex has been deleted!!!!\n");
            for (int i = 0; i < temp_n; i++) {
                degree[i] = degree_temp[i];
            }
            turn_to_oldgraph();
            continue;
        }

        addtoS(nV[q]);
        memset(ins, 0, n * sizeof(bool));
        memset(del, 0, n * sizeof(bool));
        memset(notadj, 0, n * sizeof(int));
        if (bb_s_plex_new(n, start, finish, timeflag))
        {
            if (LB < ans) ans = LB;
            if (timeflag == false)
                return 0;
        }

        for (int i = 0; i < temp_n; i++) {
            degree[i] = degree_temp[i];
        }
        turn_to_oldgraph();
    }
    return LB;
}

int fflag[MaxN] = { 0 };
int best_iterative_query(int q) {
    clock_t start, finish = 0;
    start = clock();
    bool timeflag = true;
    int y = n; S = 0;
    storeoldgraph();

    bool* tempnxt = new bool[n];
    memset(tempnxt, 0, sizeof(bool) * n);
   
    int temp_n = n;
    //这里是不是应该改成n-k?
    for (int i = 0; i < temp_n - K; ++i) {
        S++;
        cout << "S: " << S << endl;
        cout << "S + K: " << S + K << endl;
        LB = K + S;

        memset(outcore, 0, n * sizeof(bool));
        second_order_reduction(S, LB);
        seek_remain_graph();

     /*   if (S == 1)
        {
            for (int j = 0; j < nnodes; j++)
            {
                if (exist_oldid[j] == 0)
                    fflag[j] = 1;
            }
        }

        else if (S == 2)
        {
            for (int j = 0; j < nnodes; j++)
            {
                if (exist_oldid[j] && fflag[j] == 1)
                    cout << j << endl;
            }
        }
      */


        delvex.clear(); svex.clear();
        if (exist_oldid[q] == 0)
        {
            printf("Query vertex has been deleted!!!!\n");
            for (int i = 0; i < temp_n; i++) {
                degree[i] = degree_temp[i];
            }
            turn_to_oldgraph();
            continue;
        }

        memset(ins, 0, n * sizeof(bool));
        memset(del, 0, n * sizeof(bool));
        memset(notadj, 0, n * sizeof(int));

        addtoS(nV[q]);

        if (bb_s_plex_new(n, start, finish, timeflag))
        {
            if (timeflag == false)
                return 0;
            break;
        }

        for (int i = 0; i < temp_n; i++) {
            degree[i] = degree_temp[i];
        }
        turn_to_oldgraph();
    }
    return K + S;
}

int best_iterative_multi_query(vector<int> Q) {
    clock_t start, finish = 0;
    start = clock();
    bool timeflag = true;
    int y = n; S = 0;
    storeoldgraph();

    bool* tempnxt = new bool[n];
    memset(tempnxt, 0, sizeof(bool) * n);
    int temp_n = n;
    //这里是不是应该改成n-k?
    for (int i = 0; i < temp_n - K; ++i) {
        S++;
        cout << "S: " << S << endl;
        cout << "S + K: " << S + K << endl;
        LB = K + S;

        memset(outcore, 0, n * sizeof(bool));
        second_order_reduction(S, LB);
        seek_remain_graph();

        //for (int j = 0; j < nnodes;j++)
        //{
        //    if (exist_oldid[j])
        //        cout << j << endl;
        //}
        delvex.clear(); svex.clear();
        bool del_flag = false;
        for (auto q : Q)
        {
            if (!exist_oldid[q])
            {
                printf("Query vertex has been deleted!!!!\n");
                for (int i = 0; i < temp_n; i++) {
                    degree[i] = degree_temp[i];
                }
                turn_to_oldgraph();
                del_flag = true;
                break;
            }
        }
        if (del_flag) continue;

        memset(ins, 0, n * sizeof(bool));
        memset(del, 0, n * sizeof(bool));
        memset(notadj, 0, n * sizeof(int));

        for (auto q : Q) addtoS(nV[q]);

        if (bb_s_plex_new(n, start, finish, timeflag))
        {
            if (timeflag == false)
                return 0;
            break;
        }

        for (int i = 0; i < temp_n; i++) {
            degree[i] = degree_temp[i];
        }
        turn_to_oldgraph();
    }
    return K + S;
}

int best_iterative_topr(long long r) {
    clock_t start, finish = 0;
    start = clock();
    bool timeflag = true;
    int y = n; S = 0;
    long long c = 0;
    storeoldgraph();

    bool* tempnxt = new bool[n];
    memset(tempnxt, 0, sizeof(bool) * n);
    int temp_n = n;
    topr_ans.clear();
    //这里是不是应该改成n-k?
    for (int i = 0; i < temp_n - K; ++i) {
        S++;
        cout << "S: " << S << endl;
        cout << "S + K: " << S + K << endl;
        LB = K + S;

        memset(outcore, 0, n * sizeof(bool));
        memset(exist_oldid, 1, n * sizeof(int));
        second_order_reduction(S, LB);
        seek_remain_graph();

        //for (int ii = 0; ii < n; ii++)
        //{
        //    printf("degree %d: %d\n", ii, degree[ii]);
        //}

        delvex.clear(); svex.clear();


        memset(ins, 0, n * sizeof(bool));
        memset(del, 0, n * sizeof(bool));
        memset(notadj, 0, n * sizeof(int));


        if (bb_s_plex_topr(n, start, finish, timeflag, r, c) || c >= r)
        {
            if (timeflag == false)
                return 0;
            //printf("前%d个最小%d-core是：\n", r, K);
          /*  for (int i = 0; i < topr_ans.size(); i++)
            {
                for (auto j : topr_ans[i]) cout << oriID[j] << " ";
                cout << endl;
            }*/
            return S + K;
        }
        for (int i = 0; i < n; i++) {
            degree[i] = degree_temp[i];
        }
        turn_to_oldgraph();
    }
  
    return 0;
}

int best_iterative_weight(int q, double &ub) {
    clock_t start, finish = 0;
    start = clock();
    bool timeflag = true;
    int y = n; S = 0;
    long long c = 0;
    storeoldgraph();

    bool* tempnxt = new bool[n];
    memset(tempnxt, 0, sizeof(bool) * n);
    int temp_n = n;
    topr_ans.clear();
    bool flag = false;
    //这里是不是应该改成n-k?
    for (int i = 0; i < temp_n - K; ++i) {
        S++;
        cout << "S: " << S << endl;
        cout << "S + K: " << S + K << endl;
        LB = K + S;

        memset(outcore, 0, n * sizeof(bool));
        memset(exist_oldid, 1, n * sizeof(int));
        second_order_reduction(S, LB);
        seek_remain_graph();


       
        cout << endl;
        //这里要加一步来更新weight

        delvex.clear(); svex.clear();
        if (exist_oldid[q] == 0)
        {
            printf("Query vertex has been deleted!!!!\n");
            for (int i = 0; i < temp_n; i++) {
                degree[i] = degree_temp[i];
            }
            turn_to_oldgraph();
            continue;
        }
        addtoS(nV[q]);

        memset(ins, 0, n * sizeof(bool));
        memset(del, 0, n * sizeof(bool));
        memset(notadj, 0, n * sizeof(int));

        flag = false;
        bb_s_plex_weight(n, start, finish, timeflag, flag, ub, nV[q]);
        if (timeflag == false) return 0;
        if (flag) return S + K;

        for (int i = 0; i < n; i++) {
            degree[i] = degree_temp[i];
        }
        turn_to_oldgraph();
    }

    return 0;
}


void test_head(int n)
{
    cout << "test!" << endl;
    for (int i = 0; i < n; i++)
    {
        printf("%d:", i);
        for (int e = head[i]; ~e; e = nxt[e])
        {
            printf("%d ", to[e]);
        
        }
        cout << endl;
    }
}


int main(int argc, char* argv[]) {

    vector<int> kset = { stoi(argv[3]) };//3};
    ofstream outfile;

    //outfile.open("/home/zhangqifan/min_k_core/dataset/IBB_result.txt");
    //fileset = readtxt_name("/home/zhangqifan/min_k_core /dataset/name.txt");
    //fileset.clear();
    //fileset.push_back("/home/zhangqifan/min_k_core/dataset/Livejournal.txt");
    int method = 3;

    //outfile.open(argv[5]);//"/home/zhangqifan/min_k_truss/dataset/IBB_result.txt");           
    //fileset = readtxt_name("/home/heany/GX/zhangqifan/dataset_IBB/name.txt");
    fileset.clear();
    fileset.push_back(argv[1]);//"/home/zhangqifan/min_k_truss/dataset/Test.txt");

    if (method == 1) //BNB---------unweight
    {
        for (int turn = 0; turn < kset.size(); turn++)
        {
            K = kset[turn];
            for (auto it = fileset.begin(); it != fileset.end(); ++it)
            {
                refresh();
                clock_t start, finish, start1, finish1;
                string filename = *it;

                read(filename);
                //read_ratio(filename, 1);


                start = clock();
                start1 = clock();
                //oldpreprocess();
                //test_head(nnodes);
                oldpreprocess();
                //test_head(n);
                //for (int i = 0; i < nnodes; i++)
                //{
                //    if (exist_oldid[i])
                //        cout << i << endl;
                //}
                finish1 = clock();


                //nonpreprocess();
                int cnt = best_iterative();
                //vector<int> ans_temp;
                //for (int ii = 0; ii < topr_ans[0].size(); ii++)
                //{
                //    ans_temp.push_back(topr_ans[0][ii]);
                //}
                //sort(ans_temp.begin(), ans_temp.end());
                //ofstream outfile;
                //outfile.open("/home/ljw/zqf/dataset/min_k_core.txt");
                //for (int ii = 0; ii < ans_temp.size(); ii++)
                //{
                //    for (int jj = ii + 1; jj < ans_temp.size(); jj++)
                //    {
                //        if (check_neighbor(ans_temp[ii], ans_temp[jj]))
                //            outfile << ans_temp[ii] << " " << ans_temp[jj] << endl;
                //    }
                //}
                cout << "Solution4:iterative is: " << cnt << endl;
                finish = clock();
                cout << "the running time of the iterative based algorithm is:" << (double)(finish - start) / CLOCKS_PER_SEC << " s";
                outfile << "K is :" << K << " " << "and the iterative running time of" << " " << filename << "is:" << " " << (double)(finish - start) / CLOCKS_PER_SEC << "s\n";
                outfile << "The size of the min_k_core is:" << " " << cnt << endl;


            }
            outfile << endl << endl;
        }
    }
    else if (method == 2) //Enumeration---------vertex weight
    {
        string filename_weight = argv[2];// "/home/zhangqifan/min_k_truss/dataset/Test_w.txt";
        int q = stoi(argv[4]);//7;
        for (int turn = 0; turn < kset.size(); turn++)
        {
            K = kset[turn];
            for (auto it = fileset.begin(); it != fileset.end(); ++it)
            {
                refresh();
                clock_t start, finish, start1, finish1;
                string filename = *it;

                read(filename);
                read_weight(filename_weight);
                //read_ratio(filename, 1);

                start = clock();
                oldpreprocess();
                if (!exist_oldid[q])
                {
                    printf("query has been deleted!\n");
                    return 0;
                }
                q = nV[q];
             /*   for (int i = 0; i < nnodes; i++)
                {
                    if (exist_oldid[i])
                    {
                        vertex_weight[nV[i]] = vertex_weight[i];
                    }
                }*/
                svex.clear();
                memset(del, 0, sizeof(int) * MaxN);
                memset(ins, 0, sizeof(int) * MaxN);
                double ub = minw_heu(q);
                addtoS_weight_enum(q);
                svex.push_back(q);
                cout << degree[q] << endl;
                int ans_size = n;
                
                enum_weight(q, -1, ub, start, 1, ans_size);
                finish = clock();
                cout << "minimum weight is: " << ub << endl;

                outfile << "K is :" << K << " " << "and the running time of" << " " << filename << "is:" << " " << (double)(finish - start) / CLOCKS_PER_SEC << "s\n";
                outfile << "The minimum weight is:" << " " << ub << endl;
                outfile<< "The size is:" << " " << ans_size << endl;
                outfile << "min size is:" << " " << ans_size << endl;
                cout << "The size is:" << " " << ans_size<< endl;
                cout << "min size is" << " " << minsize << endl;
                cout << "the running time is:" << (double)(finish - start) / CLOCKS_PER_SEC << " s";
            }
            outfile << endl << endl;
        }
    }
    else if (method == 3)
    {
        string filename_weight = argv[2];// "/home/zhangqifan/min_k_truss/dataset/Test_w.txt";
        int q = stoi(argv[4]);//7;
        for (int turn = 0; turn < kset.size(); turn++)
        {
            K = kset[turn];
            for (auto it = fileset.begin(); it != fileset.end(); ++it)
            {
                refresh();
                clock_t start, finish, start1, finish1;
                string filename = *it;

                read(filename);
                read_weight(filename_weight);
                //read_ratio(filename, 1);

                start = clock();
                oldpreprocess();
                if (!exist_oldid[q])
                {
                    printf("query has been deleted!\n");
                    return 0;
                }
               /* for (int i = 0; i < nnodes; i++)
                {
                    if (exist_oldid[i])
                    {
                        vertex_weight[nV[i]] = vertex_weight[i];
                    }
                }*/



                double ub = MaxN;
                int cnt = best_iterative_weight(q, ub);
                


                cout << "minimum weight is: " << ub << endl;
                cout << "size is: " << cnt << endl;

                outfile << "K is :" << K << " " << "and the running time of" << " " << filename << "is:" << " " << (double)(finish - start) / CLOCKS_PER_SEC << "s\n";
                outfile << "The minimum weight is:" << " " << ub << endl;
            }
            outfile << endl << endl;
        }
    }
    outfile.close();
    return 0;
}



