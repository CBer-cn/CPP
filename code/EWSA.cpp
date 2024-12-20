#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <random>
#include <functional>
#include <ctime>
#include <random>
#include<iomanip>
#include <cstring>
#include <chrono>
#include <iostream>
#include <ratio>
#include <sys/types.h>
#include <dirent.h>
#include <iomanip>
#include <string>
#include <cfloat>
#include <sys/types.h> 
#include<algorithm>
#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <random>
#include <functional>
#include <ctime>
#include <random>
#include<iomanip>
#include <cstring>
#include <chrono>
#include <iostream>
#include <ratio>
#include <sys/types.h>
#include <dirent.h>
#include <iomanip>
#include <string>
#include <cfloat>
#include <sys/stat.h> ¡¡
#include <sys/types.h> 
#include<algorithm>
using namespace std;
clock_t start,finish;
double best_time;
int is_change=1; 
int seed=time(0)%57555;
int min_cnt=20;
//int seed=1;
mt19937 pseudoRandNumGen;
void initRand(int seed) { pseudoRandNumGen = mt19937(seed); }
auto real_rand = bind(uniform_real_distribution<double>(0,1),mt19937(seed));
vector<int> getTopKIndices(vector<int>& nums, int k) {
    vector<pair<int, int>> indexedNums(nums.size());
    for (int i = 0; i < nums.size(); ++i) {
        indexedNums[i] = make_pair(nums[i], i);
    }

    partial_sort(indexedNums.begin(), indexedNums.begin() + k, indexedNums.end(),greater<pair<int, int>>());

    vector<int> topKIndices(k);
    for (int i = 0; i < k; ++i) {
        topKIndices[i] = indexedNums[i].second;
    }
    return topKIndices;
}
void init_sol(vector<vector<int>> &benefit,vector<vector<int>> &C,vector<vector<int>> &pos,vector<int> &cluster_size,int n)
{
	for(int i=0;i<n;i++)
	{
		while(C[i].size()>0)
		{
			C[i].pop_back();
		}
	}
	for(int i=0;i<n;i++)
	{
		C[i].push_back(i);
//		C[i][0]=i;
		pos[i][0]=i;
		pos[i][1]=0;
		cluster_size[i]=1;
//		cout<<"111111111111111"<<endl;
	}
}
void Move(vector<vector<int>> &C,int c,int v,vector<vector<int>> &pos,vector<int> &cluster_size)
{
	/*delete*/
//	cluster_size[pos[v][0]]--;
//	C[pos[v][0]][pos[v][1]]=C[pos[v][0]][cluster_size[pos[v][0]]];
	pos[C[pos[v][0]][pos[v][1]]][1]=pos[v][1];
//	C[pos[v][0]].pop_back();
	/*add*/
	pos[v][0]=c;
//	pos[v][1]=cluster_size[c];
//	C[c].push_back(v);
//	cluster_size[c]++;
}


int tabu_search(int g_cnt,vector<vector<int>> &benefit,int &c,int v,int old_c,vector<vector<int>> &weight)
{
	int delta=0,delta_min=99999999;
	c=-1;
	int size;
	int now_w=weight[old_c][v];
	for(int i=0;i<g_cnt;i++)
	{
		if(old_c==i)
		{
			continue;
		}
		delta=weight[i][v];
		if(delta<delta_min)
		{
			c=i;
			delta_min=delta;
		}
	}
	delta_min-=now_w;
	return delta_min;			
}					
void updta_w(int old_c,int new_c,int v,vector<vector<int>> &weight,vector<vector<int>> &benefit,int n)
{
	int i=0;
	int limit=n-2;
	for(i=0;i<limit;i+=3)
	{
		weight[old_c][i]-=benefit[v][i];
		weight[old_c][i+1]-=benefit[v][i+1];
		weight[old_c][i+2]-=benefit[v][i+2];
	}
	for(;i<n;i++)
	{
		weight[old_c][i]-=benefit[v][i];
	}	
	for(i=0;i<limit;i+=3)
	{
		weight[new_c][i]+=benefit[v][i];
		weight[new_c][i+1]+=benefit[v][i+1];
		weight[new_c][i+2]+=benefit[v][i+2];
	}
	for(;i<n;i++)
	{
		weight[new_c][i]+=benefit[v][i];
	}
}

int find_v_in_best(vector<vector<int>> &weight,vector<int> &pos,int g_cnt,int n)
{
	int best_cnt=0;
//	cout<<"yes"<<endl;
	for(int i=0;i<n;i++)
	{
		int g_best=0,w_best=weight[0][i];
		for(int j=1;j<g_cnt;j++)
		{
			if(w_best>weight[j][i])
			{
				g_best=j;
				w_best=weight[j][i];
			}
		}
		if(g_best==pos[i])
		{
			best_cnt++;
		}
	}
	return best_cnt;
}
double ResSA(int k,int k_relax,int &k_,int &f,int &f_best,
			vector<int> &pos,vector<int> &pos_best,vector<vector<int>> &benefit,int n,vector<vector<int>> &weight,
			int T_max,int T_mid,int T_min,int per,int tag,int sw)
{
	double T,delta;
	double n_in=n*per;
	int v,c,accept; 
	int f_temp=0;
	T=T_max;
	double T_low;	
	int change_iter=0,is_done=0; 
//	while(T>T_min)
	for(int i=0;i<1000;i++)
	{
		change_iter++;
//		cout<<"T:"<<T<<endl;
		int g_cnt;
//		if(T>T_mid)
//		{
//			g_cnt=k;
//		}
//		else
//		{ 
//			g_cnt=k_relax;
//		}
		g_cnt=k_relax; 
		if(g_cnt==1)
		{
			continue;
		}
		vector<double> age(n,0);
		n_in=per*n;
		for(int j=0;j<n_in;j++)
		{
			v=pseudoRandNumGen()%n;
			accept=0;
			delta=tabu_search(g_cnt,benefit,c,v,pos[v],weight);
			if(delta<=0)
			{
				accept=1; 	 
			}
			else
			{
				double T_tmp=T*pow(1.01,age[v]*(sw^1));
				double p=min((T-T_min)/(T_max-T_min),age[v]/per)*sw;
				if(real_rand()<exp(-delta/T_tmp)||real_rand()<p)
				{
					accept=1;
				}
				
//				double p=min((T-T_min)/(T_max-T_min),age[v]/per)*sw;
//				if(real_rand()<exp(-delta/T)||real_rand()<p)
//				{
//					accept=1;
//				}
				
//				if(real_rand()<exp(-delta/T)||real_rand()<p)
//				{
//					accept=1;
//				}
				
//				double T_tmp=T*pow(1.01,age[v]);
//				if(real_rand()<exp(-delta/T_tmp))
//				{
//					accept=1;
//				}
			}
			if(accept==1)
			{	
				updta_w(pos[v],c,v,weight,benefit,n);	
				pos[v]=c;
				f+=delta;
				if(f<f_best)
				{
//					finish=clock();
//					best_time=(double)(finish-start)/CLOCKS_PER_SEC;
					pos_best=pos;
					f_best=f;
				}
				if(f<f_temp || i==1)
				{
					f_temp=f;
					change_iter=0;
					if(tag==1)
					{
						T_low=T;
					} 
				}
//				if(f<f_tmp)
//				{
//					f_tmp=f;
//					C_tmp=C;
//				}
				age[v]=0;
			}
			else
			{
				age[v]++;
			}	
		}
//		if(is_done==1)
//		{
//			cout<<"hard "<<"i:"<<i<<"   f_temp:"<<f_temp<<endl;
//			break;
//		}
//		if(change_iter>=20 )
//		{
//			cout<<"i:"<<i<<"   f_temp:"<<f_temp<<endl;
////			cout<<f_temp<<endl;
////			break;
//			T=0.001;
//			is_done=1;
//		}
		
		if(T<T_min)
		{
//			cout<<"i:"<<i<<endl;
//			change_iter=0;
			break;
		}
		T*=0.98;
	}
	cout<<"f_tmp: "<<f_temp<<endl;
	if(tag==1)
	{
		return T_low;
	}
	else
	{
		return 6.6;
	}
}


void rand_sol(int &k,int &k_relax,int n,vector<int> &pos,int &f,
				vector<vector<int>> &benefit,vector<vector<int>> &weight)
{
	int g;
	f=0;
	vector<int> cluster_size(n,0);
	vector<vector<int>> C(k_relax);
	for(int i=0;i<n;i++)
	{
		g=pseudoRandNumGen()%k;
		C[g].push_back(i);
		pos[i]=g;
		for(int j=0;j<cluster_size[g];j++)
		{
			f+=benefit[i][C[g][j]];
		}
		cluster_size[g]++;
	}
	for(int i=0;i<k_relax;i++)
	{
		weight[i].resize(n);
		for(int j=0;j<n;j++)
		{
			weight[i][j]=0;
			for(int ii=0;ii<cluster_size[i];ii++)
			{
				weight[i][j]+=benefit[j][C[i][ii]];
			}
		}
	}
}

void restart(int k,int k_relax,int n,vector<int> &pos,int &f,vector<vector<int>> &benefit)
{
	int g;
	f=0;
	vector<vector<int>> C(k_relax);
	vector<int> cluster_size(k_relax,0);
	for(int i=0;i<n;i++)
	{
		g=pseudoRandNumGen()%k;
		C[g].push_back(i);
		pos[i]=g;
		for(int j=0;j<cluster_size[g];j++)
		{
			f+=benefit[i][C[g][j]];
		}
		cluster_size[g]++;
	}
}

void GetFileNames(string path,vector<string>& filenames)
{
    DIR *pDir;
    struct dirent* ptr;
    if(!(pDir = opendir(path.c_str())))
        return;
    while((ptr = readdir(pDir))!=0) {
        if (strcmp(ptr->d_name, ".") != 0 && strcmp(ptr->d_name, "..") != 0)
            filenames.push_back(ptr->d_name);
    }
    closedir(pDir);
}

void re_sol(int &k,int &k_relax,int n,vector<int> &pos,
				vector<vector<int>> &benefit,vector<vector<int>> &weight,int &f)
{
	vector<int> k_TB(k_relax,0);
	for(int i=0;i<n;i++)
	{
		k_TB[pos[i]]++;
	}
//	cout<<"k_TB: ";
//	for(int i=0;i<k_relax;i++)
//	{
//		cout<<k_TB[i]<<" ";
//	}
//	cout<<endl;
	int g;
	vector<int> cluster_size(k_relax);
	vector<vector<int>> C(k_relax);
	for(int i=0;i<n;i++)
	{
		C[pos[i]].push_back(i);
	}
	for(int i=0;i<k_relax;i++)
	{
		cluster_size[i]=C[i].size();
	}
	for(int i=0;i<k_relax;i++)
	{
		weight[i].resize(n);
		for(int j=0;j<n;j++)
		{
			weight[i][j]=0;
			for(int ii=0;ii<cluster_size[i];ii++)
			{
				weight[i][j]+=benefit[j][C[i][ii]];
			}
		}
	}
	f=0;
//	cout<<"11111"<<endl;
	for(int i=0;i<n;i++)
	{
		f+=weight[pos[i]][i];
	}
	cout<<f<<endl;
	f/=2;
}		

int TSA(int g_cnt,vector<int> &pos,vector<vector<int>> &benefit,int n,vector<vector<int>> &weight,
			int T_max,int per)
{
	double T,delta;
	double n_in=n*per;
	int v,c,accept;
	T=T_max;
	vector<int> age(n,0);
	for(int j=0;j<n_in;j++)
	{
		accept=0;
		v=pseudoRandNumGen()%n;
		delta=tabu_search(g_cnt,benefit,c,v,pos[v],weight);
		if(delta<=0)
		{
			accept=1;
		}
		else
		{
			double T_tmp=T*pow(1.01,age[v]);
			if(real_rand()<exp(-delta/T_tmp))
			{
				accept=1;
			}
			
			
		}
		if(accept==1)
		{
			updta_w(pos[v],c,v,weight,benefit,n);
			pos[v]=c;
			age[v]=0;
		}
		else
		{
			age[v]++;
		}	
	}
	return find_v_in_best(weight,pos,g_cnt,n);
}

void InitTemperature(double &T_max,double &T_mid,double &T_min,int g_cnt,vector<vector<int>> &benefit,int n)
{
	int obj1=n*0.8;
	int obj2=n*0.9; 
	int obj3=n*1; 
	double T_l=2000,T_r=0;
	
	while(1)
	{
		vector<vector<int>> weight(g_cnt);
		vector<int> pos(n,0);
		int f;
		rand_sol(g_cnt,g_cnt,n,pos,f,benefit,weight);
		int obj=TSA(g_cnt,pos,benefit,n,weight,T_max,100);
//		cout<<obj<<endl;
		if(obj-obj1<5 &&obj-obj1>-5)
		{
			break;
		}
		else if(obj-obj1>=5)
		{
			T_r=T_max;
			T_max=(T_l+T_r)/2;
		}
		else if(obj-obj1<=-5)
		{
			T_l=T_max;
			T_max=(T_l+T_r)/2;
		}
//		cout<<T_l<<"  "<<T_r<<endl;
		if(abs(T_l-T_r)<1)
		{
			break;
		}
	}
	
	T_l=2000,T_r=0;
	while(1)
	{		
		vector<vector<int>> weight(g_cnt);
		vector<int> pos(n,0);
		int f;
		rand_sol(g_cnt,g_cnt,n,pos,f,benefit,weight);
		int obj=TSA(g_cnt,pos,benefit,n,weight,T_mid,100);
//		cout<<obj<<endl;
		if(obj-obj2<5 &&obj-obj2>-5)
		{
			break;
		}
		else if(obj-obj2>=5)
		{
			T_r=T_mid;
			T_mid=(T_l+T_r)/2;
		}
		else if(obj-obj2<=-5)
		{
			T_l=T_mid;
			T_mid=(T_l+T_r)/2;
		}
//		cout<<T_l<<"  "<<T_r<<endl;
		if(abs(T_l-T_r)<1)
		{
			break;
		}
	}
	
	int f;
	int g_cnt_;
	vector<vector<int>> weight(g_cnt);
	vector<int> pos(n,0);
	rand_sol(g_cnt,g_cnt,n,pos,f,benefit,weight);
	vector<int> pos_best=pos;
	int f_best=f;
	
//	double T_min1=ResSA(g_cnt,g_cnt,g_cnt_,f,f_best,pos,pos_best,benefit,n,weight,T_max,T_mid,T_min,100,1,0);
	double T_min1=ResSA(g_cnt,g_cnt,g_cnt_,f,f_best,pos,pos_best,benefit,n,weight,T_mid,T_mid,T_min,100,1,0);
	T_min=T_min1;
	if(T_min>=10)
	{
		int tcnt=T_min/10;
		T_min=tcnt*10;
	}
	else if(T_min<10 &&T_min>=5)
	{
		T_min=5;
	}
	else if(T_min<5&&T_min>=2)
	{
		T_min/=2;
	}
	else if(T_min<2&&T_min>=1)
	{
		T_min=1;
	}
	else if(T_min<1&&T_min>=0.5)
	{
		T_min=0.5;
	}
	else
	{
		T_min=0.1;
	}
}

void ReduceGroup(vector<int> &pos,int k,int n)
{
	int k_now=-1;
	for(int i=0;i<n;i++)
	{
		if(pos[i]>k_now)
		{
			k_now=pos[i];
		}
	}
	k_now=k_now+1;
	vector<int> k_TB(k_now,0);
	for(int i=0;i<n;i++)
	{
		k_TB[pos[i]]++;
	}
	vector<vector<int>> clusters(k_now);
	for(int i=0;i<n;i++)
	{
		clusters[pos[i]].push_back(i);
	}
	vector<int> topKIndicesN = getTopKIndices(k_TB, k_TB.size());
	int cur=0;
	int save=0;
	for(int i=0;i<clusters.size();i++)
	{
		if(cur<k)
		{
			save+=clusters[topKIndicesN[i]].size();
			for(int j=0;j<clusters[topKIndicesN[i]].size();j++)
			{
				pos[clusters[topKIndicesN[i]][j]]=cur;
			}
			cur++;
		}
		else
		{
			for(int j=0;j<clusters[topKIndicesN[i]].size();j++)
			{
				pos[clusters[topKIndicesN[i]][j]]=pseudoRandNumGen()%k;
			}
		}
	}
//	cout<<"ReduceGroup save:"<< save<<endl;
}
void Perturb(vector<int> &pos,double p,int k,int n)
{
	for(int i=0;i<n;i++)
	{
		if(real_rand()>=p)
		{
			pos[i]= pseudoRandNumGen()%k;
		}
	}
}
void Oscillate(int k,int k_relax,double &T_init,double T_max,double T_mid,vector<int> &pos,vector<int> &pos_best,
			vector<vector<int>> &benefit,int n,vector<vector<int>> &weight,int &f,int f_best,int &I_in,int &sw)
{
	
//	rand_sol(vector<vector<int>> &C,int g_cnt,int n,vector<int> &cluster_size,vector<vector<int>> &pos,double &f,
//				vector<vector<int>> &benefit,vector<vector<int>> &weight)
	if(sw==0)
	{
//		cout<<"block!"<<endl;
		T_init=T_mid;
//		sw=0;			
//		cout<<"fine-grained "<<endl;
		cout<<"fine-grained k£º"<<k<<endl;
		f=f_best;
//		vector<int> k_TB(k_relax,0);
//		for(int i=0;i<n;i++)
//		{
////			cout<<pos[i]<<endl;
//			k_TB[pos_best[i]]++;
//		}
//		cout<<"k_TB: ";
//		for(int i=0;i<k_relax;i++)
//		{
//			cout<<k_TB[i]<<" ";
//		}
//		cout<<endl;
		pos=pos_best;
//		Perturb(pos,0.9,k,n);
		re_sol(k,k_relax,n,pos,benefit,weight,f);
//		I_in=50;
		I_in=100;
//		if(n/k/2>100)
//		{
//			I_in=n/k/2;
//		}
//		else
//		{
//			I_in=100;
//		}
	}
	else
	{
		T_init=T_max;
		cout<<"coarse-grained k£º"<<k<<endl;
		ReduceGroup(pos,k,n);
//		Perturb(pos,0.8,k,n);
		re_sol(k,k_relax,n,pos,benefit,weight,f);
//		rand_sol(k,k_relax,n,pos,f,benefit,weight);
		I_in=100;
//		if(n/k/2>100)
//		{
//			I_in=n/k/2;
//		}
//		else
//		{
//			I_in=100;
//		}
	}
//	I_in=100;
}


double initTime(int n)
{
//	cout<<"n:"<<n<<endl;
	if(n<=300)
	{
		return 200;
	}
	else if(300<n&&n<=500)
	{
		return 500;
	}
	else if(500<n&&n<=800)
	{
		return 1000;
	}
	else if(800<n&&n<=1000)
	{
		return 2000;
	}
	else if(1000<n&& n<=1500)
	{
		return 4000;
	}
	else if(1500<n&&n<=2500)
	{
		return 10000;
	}
	else if(2500<n)
	{
		return 20000;
	}
	return 0;
}

int Caculate_delta(vector<int> &Ci,vector<vector<int>> &benefit,int v)
{
	int delta=0;
	for(int i=0;i<Ci.size();i++)
	{
		delta+=benefit[v][Ci[i]];
	}
	return delta;
}

int GreedyInit(vector<vector<int>> &benefit,int n) 
{
	int k=1,id;
	vector<vector<int>> C(n);
	vector<int> tabu_list(n,0);
	for(int i=0;i<n;i++)
	{
		int v=pseudoRandNumGen()%n;
		while(tabu_list[v]==1)
		{
			v=pseudoRandNumGen()%n;
		}
		tabu_list[v]=1;
		int delta_max=999999;
		id=-1;
		for(int j=0;j<k;j++)
		{
			int delta=Caculate_delta(C[j],benefit,v);
			if(delta_max>=delta)
			{
				delta_max=delta;
				id=j;
			}
		}
		C[id].push_back(v);
		if(id==k-1)
		{
			k++;
		}
	}
	return k;
}
void Updatek(int &k_best,int &k_min,int &k_)
{
	k_best=k_;
//	if(k_<k_min)
//	{
//		k_min=k_;
//	}
}

int AdjustPos(int &k_relax,vector<int> &pos,int &n)
{
	vector<int> k_TB(k_relax,0);
	for(int i=0;i<n;i++)
	{
		k_TB[pos[i]]++;
	}
	int cur=0;
	for(int i=0;i<k_relax;i++)
	{
		if(k_TB[i]>0)
		{
			for(int j=0;j<n;j++)
			{
				if(pos[j]==i)
				{
					pos[j]=cur;
				}
			}
			cur++;
		}
	}
	return cur;
} 
vector<vector<int>> Generate_map(vector<int> &pos,int &n,int k_relax)
{
	vector<vector<int>> p(k_relax);
	for(int i=0;i<n;i++)
	{
		p[pos[i]].push_back(i);
	}
	vector<vector<int>> map(n,vector<int>(n, 0));
	int size=p.size();
//	cout<<"size: ";
	for(int i=0;i<size;i++)
	{
		int size_=p[i].size();
//		cout<<size_<<" ";
		for(int j=0;j<size_;j++)
		{
			for(int k=j;k<size_;k++)
			{
				map[p[i][j]][p[i][k]]=1;
				map[p[i][k]][p[i][j]]=1;
			}
		}
	}
//	cout<<endl;
	return map;
}

vector<vector<int>> Generate_cross_matrix(vector<int> &sequences,vector<vector<vector<int>>> &maps,int &n)
{
	vector<vector<int>> cross(n,vector<int>(n, 1));
	int size=sequences.size();
	for(int i=0;i<size;i++)
	{
		int id=sequences[i];
		for(int j=0;j<n;j++)
		{
			for(int k=0;k<n;k++)
			{
				cross[j][k]=cross[j][k]&maps[id][j][k];
			}
		}
	}
	return cross;
}

vector<vector<int>> Generate_small_clusters(vector<vector<int>> &cross,int n)
{
	vector<vector<int>> clusters;
	vector<int> is_exist(n,0);
	for(int i=0;i<n;i++)
	{
		vector<int> cluster;
		if(is_exist[i]==1)
		{
			continue;
		}
		for(int j=0;j<n;j++)
		{
			if(cross[i][j]==1)
			{
				cluster.push_back(j);
				is_exist[j]=1;
			}
		}
		clusters.push_back(cluster);
	}
	return clusters;
} 
int Generate_new_benefit(vector<vector<int>> &benefit,vector<int> &A,vector<int> &B)
{	
	int a_size=A.size(),b_size=B.size(); 
	int sum=0;
	for(int i=0;i<a_size;i++)
	{
		for(int j=0;j<b_size;j++)
		{
			sum+=benefit[A[i]][B[j]];
		}
	}
	return sum;
}



void CoverGD(vector<int> &pos,vector<vector<int>> &clusters,vector<int> &topKIndicesN,
							int k_best,int k_max,int n,vector<int> &tabu_clusters,double &save)
{//Ëæ»úÈÅ¶¯Ê£Óàµã 
	int cur=0;
	vector<int> tabu(k_max,0);
	for(int i=0;i<clusters.size();i++)
	{
		if(tabu[pos[clusters[topKIndicesN[i]][0]]]==0 &&cur<k_best)
		{
			tabu_clusters[topKIndicesN[i]]=1;
			save+=clusters[topKIndicesN[i]].size();
			cur++;
			tabu[pos[clusters[topKIndicesN[i]][0]]]=1;
		}
	}
//	cout<<"save vertex: "<<save<<endl;
}

int population_hybridization(vector<int> &pos1,vector<int> &pos2,int &n,
								vector<vector<int>> &benefit,int &k_best,int k_max)
{
	vector<vector<vector<int>>> maps(2);
	maps[0]=Generate_map(pos1,n,k_max);
//	cout<<"------ map0 over ------"<<endl;
	maps[1]=Generate_map(pos2,n,k_max);
//	cout<<"------ map1 over ------"<<endl;
//	cout<<"------ Generate_map over ------"<<endl;
	vector<int> sequences;
	sequences.push_back(0);
	sequences.push_back(1);

	vector<vector<int>> cross=Generate_cross_matrix(sequences,maps,n);
//	cout<<"------ cross over ------"<<endl;
	
	vector<vector<int>> clusters=Generate_small_clusters(cross,n);
	
//	cout<<"------ slusters over ------"<<endl;
	
	vector<int> cw(clusters.size());
	
	vector<int> cn(clusters.size());
	
	vector<int> tabu(clusters.size(),0);
	
	
	for(int i=0;i<clusters.size();i++)
	{
		cn[i]=clusters[i].size();
	}

	cout<<"cn:"<<cn.size()<<endl;
	
	vector<int> tabu_clusters(cn.size(),0);
	vector<int> topKIndicesN = getTopKIndices(cn, cn.size());
	double save1=0,save2=0; 
	CoverGD(pos1,clusters,topKIndicesN,k_best,k_max,n,tabu_clusters,save1);
	CoverGD(pos2,clusters,topKIndicesN,k_best,k_max,n,tabu_clusters,save2);
	double repeatability=(save1+save2)/2/n;
	cout<<"repeatability:"<<repeatability<<endl;
	if(repeatability>0.9)
	{
		return 1;
	}
	else
	{
		return 0;
	}
}


int main(int argc, char **argv)
{
	char oppath[20];
	vector<string> files; 
	string inputpath=argv[1];//Folder where the examples are located 
	string outputpath=argv[2];//Solution file output location
//	string inputpath="data";
//	string outputpath="solution";
//	string save_time="1000";
//	int limit_time=stod(save_time);
	strcpy(oppath,outputpath.c_str());
	GetFileNames(inputpath,files);
	for(int i=0;i<files.size();i++)
	{
		cout<<files[i]<<endl;
	}
//	int isCreate = mkdir(oppath,S_IRUSR | S_IWUSR | S_IXUSR | S_IRWXG | S_IRWXO);	
	for(int file_id=0;file_id<files.size();file_id++)
	{
		for(int run_id=0;run_id<1; run_id++)
		{
			ifstream in(inputpath+"/"+files[file_id]);
			ofstream fout(outputpath+"/"+files[file_id], ios::app); 
			initRand(seed);
			
			int n;
			in>>n;	
			double limit_time=initTime(n);
			vector<vector<int>> benefit(n);	
			vector<int> pos_best(n);//solution
			int f_best=0,f;
			vector<int> pos(n);
			vector<int> cluster_size(n);
			for(int i=0;i<n;i++)
			{
				benefit[i].resize(n);
			}
			/*The second type of instance reads in*/
			for(int i=0;i<n;i++)
			{
				for(int j=i;j<n;j++)
				{
					int x;
					in>>x;
					benefit[i][j]=x;
					benefit[j][i]=x;
				}
			}	
			/*The first type of instance reads in*/
//			for(int i=0;i<n-1;i++)
//			{
//				for(int j=i+1;j<n;j++)
//				{
//					int x;
//					in>>x;
//					benefit[i][j]=x*-1;
//					benefit[j][i]=x*-1;
//				}
//			}	 
			 
			int k_=GreedyInit(benefit,n);
			k_=k_+sqrt(k_);
//			k_=20;
			int k_best=k_,k_min=k_;
			cout<<"init k:"<<k_<<endl; 
			double totaltime;
			start=clock();
			vector<int> v_list;
			double T_max=1000,T_mid=1000,T_min=0.01;
			InitTemperature(T_max,T_mid,T_min,k_best,benefit,n);
			totaltime=(double)(finish-start)/CLOCKS_PER_SEC;
			cout<<"totaltime: "<<totaltime<<endl;
//			T_min=0.1;
//			T_max=1.46,T_mid=0.5,T_min=0.1;
			cout<<"T_max:"<<T_max<<"       T_mid:"<<T_mid<<"      T_min:"<<T_min<<endl;
//			return 0;
//			continue;
			double beta=0.5;
			int flag=0;
			int I_in;
			vector<vector<int>> weight(k_best);
			rand_sol(k_best,k_best,n,pos_best,f_best,benefit,weight);
			pos=pos_best;
//			rand_sol(int &k,int &k_relax,int n,vector<int> &pos,int &f,
//				vector<vector<int>> &benefit,vector<vector<int>> &weight)
//			rand_sol(C,k,k_relax,n,cluster_size,pos,f,benefit,weight);
			int is_start;
			double T_init;
//			vector<vector<int>> C_tmp=C_best;
//			int f_tmp=f_best;
			int unimprove_time=0;
			vector<int> pos_gl=pos_best;//solution
			int f_gl=f_best,f_worse=-999999999;
			int f_best_last=f_best;
			int turn=0,k_now;
			vector<vector<int>> past_best;
			double op; 
//			vector<int> pos(n);
			while(1)
			{
				cout<<"================= "<<turn<<" =================="<<endl;
				int k,k_relax; 
				int sw;
				int maxspan=f_worse-f_gl;
				int span=f_worse-f_best;
				totaltime=(double)(finish-start)/CLOCKS_PER_SEC;
				if(maxspan==0)
				{
					op=(double)min(1.0,2*totaltime/limit_time);
				}
				else
				{
					op=(double)min(1.0,2*totaltime/limit_time)*(double)span/maxspan;
				}
				cout<<"!!!! op:"<<op<<" !!!!      span:"<<span<<"   maxspan:"<<maxspan<<endl;
				if(real_rand()<op)
				{
					sw=0;
					k=k_best+1;
					k_relax=k_best+1; 
					unimprove_time++;
				}
				else
				{
					sw=1;
					k=k_best+1;
					k_relax=k_best+1; 
				}
				
//				vector<int> pos(n);
				vector<vector<int>> weight(k_relax);
				double T_init;
				Oscillate(k,k_relax,T_init,T_max,T_mid,pos,pos_best,benefit,n,weight,f,f_best,I_in,sw);
				cout<<"f start:"<<f<<endl;
//				ResSA(k,k_relax,k_,f,f_best,pos,cluster_size,benefit,n,weight,T_init,T_mid,T_min,0,I_in,3,1,sw);
				ResSA(k,k_relax,k_,f,f_best,pos,pos_best,benefit,n,weight,T_init,T_mid,T_min,I_in,1,sw);
				cout<<"f_best:"<<f_best*-1<<endl;
//				cout<<"befor updata time:"<<(double)(finish-start)/CLOCKS_PER_SEC-totaltime<<endl;
				if(f_best>f_worse)
				{
					f_worse=f_best;
				}
				int pass=0;
				
				
				
				if(f_best_last>f_best)
				{
					
					ReduceGroup(pos_best,k_relax,n);
					for(int i=0;i<past_best.size();i++)
					{
						ReduceGroup(past_best[i],k_relax,n);
						pass=population_hybridization(pos_best,past_best[i],n,benefit,k_best,k_relax);
						if(pass==1)
						{
							break;
						}
					}
//					if(pass!=1)
//					{
//						past_best.push_back(pos_best);
//					}
					unimprove_time=0;
					f_best_last=f_best;
					k_now=AdjustPos(k_relax,pos_best,n);	
				}
				if(f_gl>f_best)
				{
					f_gl=f_best;
					pos_gl=pos_best;
					unimprove_time=-10;
					Updatek(k_best,k_min,k_now);
					finish=clock();
					best_time=(double)(finish-start)/CLOCKS_PER_SEC;
				}
				cout<<"past_best size:"<<past_best.size()<<"   pass:"<<pass<<"    unimprove_time:"<<unimprove_time<<endl;
				if(unimprove_time>15 || pass==1)
				{
					unimprove_time=0;
					if(pass!=1)
					{
						past_best.push_back(pos_best);
					} 
					
					restart(k_best,k_best,n,pos_best,f_best,benefit);
					f_best_last=f_best;
				}

				cout<<"f_gl:"<<f_gl*-1<<endl;
				cout<<"unimprove_timet:"<<unimprove_time<<endl;
				finish=clock();
				cout<<"each time:"<<(double)(finish-start)/CLOCKS_PER_SEC-totaltime<<endl;
				totaltime=(double)(finish-start)/CLOCKS_PER_SEC;
				cout<<"totaltime: "<<totaltime<<endl;
				if(limit_time<totaltime)
				{
					break;
				}
//				if(20<totaltime)
//				{
//					break;
//				}
				turn++;
			}
			cout<<"***********************************************************************************************"<<endl;
			vector<vector<int>> C_gl(k_best+1);
			for(int i=0;i<n;i++)
			{
//				fout<<"---------------------------------------"<<endl;
				fout<<pos_gl[i]<<" ";
			}
			fout<<endl;
			fout<<"---------------------------------------"<<endl;
			fout<<"f_best "<<f_gl*-1<<endl;
			fout<<"best_time "<<best_time<<endl;
			fout<<"==========================================================="<<endl;
			cout<<"f_best:"<<f_gl*-1<<endl;
		}
	}
	return 0; 
}
