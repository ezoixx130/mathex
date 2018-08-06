#include<cmath>
#include<complex>
#include<cstdio>
#include<cstring>
#include<iostream>
#include<map>
#include<vector>

using namespace std;

const long long PerfectNumber[]={0ll,6ll,28ll,496ll,8128ll,33550336ll,8589869056ll,137438691328ll,2305843008139952128ll};
const double C_Pi=acos(-1.);

long long GetRandomInt(long long x,long long y){
	long long d=y-x+1,R=rand()*32768ll*32768ll%d*32768ll%d+rand()*32768ll*32768ll%d+rand()*32768ll%d+rand()%d;
	R%=d;
	return R+x;
}

template<class Type> Type Sum(vector<Type> V){
	auto res=0ll;
	for(auto i:V)res+=i;
	return res;
}

template<class Type> Type Product(vector<Type> V){
	auto res=1ll;
	for(auto i:V)res*=i;
	return res;
}

namespace Natural{
	long long Gcd(long long x,long long y){
		return x?Gcd(y%x,x):y;
	}
	long long Lcm(long long x,long long y){
		return x/Gcd(x,y)*y;
	}
	long long Power(long long a,long long x){
		long long res=1ll;
		while(x){
			if(x&1ll)res*=a;
			a*=a;
			x>>=1;
		}
		return res;
	}
	long long PowerMod(long long a,long long x,const long long m){
		long long res=1ll;
		a%=m;
		while(x){
			if(x&1ll)(res*=a)%=m;
			(a*=a)%=m;
			x>>=1;
		}
		return res;
	}
	long long Factorial(long long x){
		auto res=1ll;
		for(auto i=2ll;i<=x;++i)res*=i;
		return res;
	}
	vector<long long> GetPrime(long long N){
		map<long long,bool> check;
		vector<long long> prime;
		for(auto i=2ll;i<=N;++i){
			if(!check[i])prime.push_back(i);
			for(long long j:prime){
				if(i*j>N)break;
				check[i*j]=true;
				if(i%j==0)break;
			}
		}
		return prime;
	}
	vector<pair<long long,long long>> PrimeFactorize(long long N){
		vector<pair<long long,long long>> res;
		auto prime=GetPrime((long long)ceil(sqrt(N)));
		for(auto i:prime){
			auto cnt=0ll;
			while(N%i==0ll){
				N/=i;
				++cnt;
			}
			if(cnt)res.push_back(make_pair(i,cnt));
		}
		if(N>1ll)res.push_back(make_pair(N,1ll));
		return res;
	}
	bool IsPrime(long long N,long long CheckTimes=80ll){
		if(N==2ll||N==3ll)return true;
		if(N==1ll||N%2ll==0ll)return false;
		auto R=N-1ll,S=0ll;
		while(R%2ll==0ll){
			R>>=1;
			++S;
		}
		for(int i=0;i<CheckTimes;++i){
			auto A=GetRandomInt(2ll,N-1ll);
			auto M=PowerMod(A,R,N);
			if(M==1ll||M==N-1ll)continue;
			auto j=0ll;
			while(j<S&&M!=N-1ll){
				M=M*M%N;
				if(M==1ll)return false;
				++j;
			}
			if(M!=N-1ll)return false;
		}
		return true;
	}
	bool IsSquare(long long N){
		long double Root=sqrt(N);
		long long Rt1=(long long)floor(Root),Rt2=(long long)ceil(Root);
		return Rt1*Rt1==N&&Rt2*Rt2==N;
	}
	namespace NumberTheoryFunction{
		long long pi(long long N){
			auto prime=GetPrime(N);
			return prime.size();
		}
		long long d(long long N){
			long long count=1ll;
			auto factor=PrimeFactorize(N);
			for(auto j:factor)count*=j.second+1ll;
			return count;
		}
		long long sigma(long long N){
			long long sum=1ll;
			auto factor=PrimeFactorize(N);
			for(auto j:factor)sum*=(Power(j.first,j.second+1ll)-1ll)/(j.first-1ll);
			return sum;
		}
		long long mu(long long N){
			if(N==1ll)return 1ll;
			auto res=1ll;
			auto factor=PrimeFactorize(N);
			for(auto j:factor)if(j.second>1ll)return 0ll;else res=-res;
			return res;
		}
		long long phi(long long N){
			auto factor=PrimeFactorize(N);
			for(auto j:factor)N=N/j.first*(j.first-1ll);
			return N;
		}
	}
}

namespace Modulus{
	long long Power(long long a,long long x,const long long m){
		long long res=1ll;
		a%=m;
		while(x){
			if(x&1ll)(res*=a)%=m;
			(a*=a)%=m;
			x>>=1;
		}
		return res;
	}
	long long Factorial(long long x,const long long m){
		auto res=1ll;
		for(auto i=2ll;i<=x;++i)(res*=i)%=m;
		return res;
	}
}

template<class Type> struct Polynomial{
	private:
		vector<Type> __element;
		bool __range_check;
		int out_of_range(){
			puts("Calling for element out of range!");
			return 0x56417ECB;
		}
		void fast_fourier_transform(int __len,vector<int> __rev,int __check){
			complex<long double> __f[__len];
			for(int i=0;i<__len;++i)__f[i]=__element[i];
			for(int i=0;i<__len;++i)if(i<__rev[i])swap(__f[i],__f[__rev[i]]);
			for(int i=1;i<__len;i<<=1){
				complex<long double> __w_n(cos(C_Pi/(long double)i),__check*sin(C_Pi/(long double)i));
				for(int j=0;j<__len;j+=(i<<1)){
					complex<long double> __w(1.,0.);
					for(int k=0;k<i;++k){
						complex<long double> __x=__f[j+k],__y=__w*__f[i+j+k];
						__f[j+k]=__x+__y;
						__f[i+j+k]=__x-__y;
						__w*=__w_n;
					}
				}
			}
			if(__check==-1)for(int i=0;i<__len;++i)__f[i]/=(long double)__len;
			for(int i=0;i<__len;++i)__element[i]=(Type)__f[i];
		}
	public:
		Polynomial(int Siz,const bool RgCk=false){
			__range_check=RgCk;
			for(int i=0;i<Siz;++i)__element.push_back(Type());
		}
		Polynomial(int Siz,Type Val,const bool RgCk=false){
			__range_check=RgCk;
			for(int i=0;i<Siz;++i)__element.push_back(Val);
		}
		Polynomial(vector<Type> &Vec,const bool RgCk=false){
			__range_check=RgCk;
			for(int i=0;i<Vec.size();++i)__element.push_back(Vec[i]);
		}
		Polynomial(Polynomial<Type> &F,const bool RgCk=false){
			__range_check=RgCk;
			for(int i=0;i<F.Size();++i)__element.push_back(F[i]);
		}
		Polynomial& operator=(vector<Type> &Vec){
			__element=Vec;
			return *this;
		}
		friend bool operator==(Polynomial<Type> &F,Polynomial<Type> &G){
			F.ShrinkToFit();
			G.ShrinkToFit();
			if(F.Size()!=G.Size())return false;
			for(int i=0;i<F.Size();++i)if(F[i]!=G[i])return false;
			return true;
		}
		friend Polynomial<Type> operator+(Polynomial<Type> &F,Polynomial<Type> &G){
			Polynomial<Type> H(max(F.Size(),G.Size()));
			for(int i=0;i<F.Size();++i)H[i]+=F[i];
			for(int i=0;i<G.Size();++i)H[i]+=G[i];
			return H;
		}
		friend Polynomial<Type> operator-(Polynomial<Type> &F,Polynomial<Type> &G){
			Polynomial<Type> H(max(F.Size(),G.Size()));
			for(int i=0;i<F.Size();++i)H[i]-=F[i];
			for(int i=0;i<G.Size();++i)H[i]-=G[i];
			H.ShrinkToFit();
			return H;
		}
		friend Polynomial<Type> operator*(Polynomial<Type> &F,Type &lambda){
			Polynomial<Type> H(F);
			for(int i=0;i<F.Size();++i)H[i]*=lambda;
			return H;
		}
		friend Polynomial<Type> operator*(Polynomial<Type> &F,Polynomial<Type> &G){
			int __size=F.Size()+G.Size()-1,__len=1,__bit=0;
			while(__len<=__size){
				__len<<=1;
				++__bit;
			}
			vector<int> __rev(__len);
			for(int i=0;i<__len;++i)__rev[i]=((__rev[i>>1]>>1)|((i&1)<<__bit-1));
			Polynomial<Type> F1=F,F2=G,H(__size);
			F1.Adjust(__size);
			F2.Adjust(__size);
			F1.fast_fourier_transform(__size,__rev,1);
			F2.fast_fourier_transform(__size,__rev,1);
			for(int i=0;i<__len;++i)H[i]=F1[i]*F2[i];
			H.fast_fourier_transform(__size,__rev,-1);
			return H;
		}
		friend Polynomial<Type> operator/(Polynomial<Type> &F,Type &lambda){
			Polynomial<Type> H(F);
			for(int i=0;i<F.Size();++i)H[i]/=lambda;
			return H;
		}
		Polynomial& operator+=(Polynomial<Type> &F){
			if(Size()<F.Size())Adjust(F.Size());
			for(int i=0;i<F.Size();++i)__element[i]+=F[i];
			return *this;
		}
		Polynomial& operator-=(Polynomial<Type> &F){
			if(Size()<F.Size())Adjust(F.Size());
			for(int i=0;i<F.Size();++i)__element[i]-=F[i];
			return *this;
		}
		Type& operator[](int x){
			if(x>=__element.size()){
				if(__range_check)throw out_of_range();else{
					while(__element.size()<=x)__element.push_back(Type());
					return __element.back();
				}
			}
			return __element[x];
		}
		void Adjust(int Siz){
			Type Emp=Type();
			while(__element.size()<Siz)__element.push_back(Emp);
			while(__element.size()>Siz)__element.pop_back();
		}
		void Print(){
			for(int i=0;i<__element.size();++i)cout<<__element[i]<<ends;
			puts("");
		}
		void Scan(int Siz=0){
			if(Siz==0)Siz=Size();
			__element.clear();
			Adjust(Siz);
			for(int i=0;i<Siz;++i)cin>>__element[i];
		}
		void SetRangeCheck(bool Sta){
			__range_check=Sta;
		}
		void ShrinkToFit(){
			Type Emp=Type();
			while(!__element.empty()&&__element.back()==Emp)__element.pop_back();
		}
		int Size(){
			return __element.size();
		}
};
