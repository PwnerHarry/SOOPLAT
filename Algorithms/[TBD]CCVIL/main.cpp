#include "Header.h"


int main(int argc, char **argv){

	RunParameter* runParam;
	unsigned funcIndex;
	Benchmarks* fp=NULL;
	
	if (argc==2){
//		specify the file of configuration
		runParam = new RunParameter(argv[1]);
	} else{
		// if no file specified, use the default file "configure.ini"
		runParam = new RunParameter();
	}
	vector<unsigned> funcToRun=runParam->functionToRun;

	
	CCVIL* algo = new CCVIL(runParam);

	for (funcIndex = 0; funcIndex < runParam->functionToRun.size(); funcIndex++ ){
		fp = generateFuncObj(runParam, funcToRun[funcIndex]);
		algo->setObjFunc(fp);
		algo->run();
		delete fp;
	}

	delete algo;
	delete runParam;

	return 0;
}


// create new object of class according to the configure.ini file
Benchmarks* generateFuncObj(RunParameter* runParam, int funcID){
	Benchmarks *fp;
	// run each of specified function in "configure.ini"
	if (funcID==1){
		fp = new F1(runParam);
	}else if (funcID==2){
		fp = new F2(runParam);
	}else if (funcID==3){
		fp = new F3(runParam);
	}else if (funcID==4){
		fp = new F4(runParam);
	}else if (funcID==5){
		fp = new F5(runParam);
	}else if (funcID==6){
		fp = new F6(runParam);
	}else if (funcID==7){
		fp = new F7(runParam);
	}else if (funcID==8){
		fp = new F8(runParam);
	}else if (funcID==9){
		fp = new F9(runParam);
	}else if (funcID==10){
		fp = new F10(runParam);
	}else if (funcID==11){
		fp = new F11(runParam);
	}else if (funcID==12){
		fp = new F12(runParam);
	}else if (funcID==13){
		fp = new F13(runParam);
	}else if (funcID==14){
		fp = new F14(runParam);
	}else if (funcID==15){
		fp = new F15(runParam);
	}else if (funcID==16){
		fp = new F16(runParam);
	}else if (funcID==17){
		fp = new F17(runParam);
	}else if (funcID==18){
		fp = new F18(runParam);
	}else if (funcID==19){
		fp = new F19(runParam);
	}else if (funcID==20){
		fp = new F20(runParam);
	}else{
		cerr<<"Fail to locate Specified Function Index"<<endl;
		exit(-1);
	}
	return fp;
}

// create new object of class with default setting
Benchmarks* generateFuncObj(int funcID){
	Benchmarks *fp;
	// run each of specified function in "configure.ini"
	if (funcID==1){
		fp = new F1();
	}else if (funcID==2){
		fp = new F2();
	}else if (funcID==3){
		fp = new F3();
	}else if (funcID==4){
		fp = new F4();
	}else if (funcID==5){
		fp = new F5();
	}else if (funcID==6){
		fp = new F6();
	}else if (funcID==7){
		fp = new F7();
	}else if (funcID==8){
		fp = new F8();
	}else if (funcID==9){
		fp = new F9();
	}else if (funcID==10){
		fp = new F10();
	}else if (funcID==11){
		fp = new F11();
	}else if (funcID==12){
		fp = new F12();
	}else if (funcID==13){
		fp = new F13();
	}else if (funcID==14){
		fp = new F14();
	}else if (funcID==15){
		fp = new F15();
	}else if (funcID==16){
		fp = new F16();
	}else if (funcID==17){
		fp = new F17();
	}else if (funcID==18){
		fp = new F18();
	}else if (funcID==19){
		fp = new F19();
	}else if (funcID==20){
		fp = new F20();
	}else{
		cerr<<"Fail to locate Specified Function Index"<<endl;
		exit(-1);
	}
	return fp;
}
