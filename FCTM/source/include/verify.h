#ifndef VERIFY_H
#define VERIFY_H

void getPD(dtype *s, int n, int *PD){
    int st=-1;
    printf("PD: ");
    dtype *stv = (dtype *) malloc(n*sizeof(dtype));
    int *sti = (int *) malloc(n*sizeof(int));
    for(int i=0;i<n;i++){
        dtype value;
        int index;
        while(st!=-1){
            value=stv[st];
            index=sti[st];
            if(value<=s[i]) break;
            st--;
        }
        if(st==-1) PD[i]=0;
        else PD[i]=i-index;
        printf("%d ", PD[i]);
        stv[++st]=s[i];
        sti[st]=i;
    }
    printf("\n");
    free(stv);
    free(sti);
    return;
}


int verification(dtype *x, int m, dtype *y, int n, int *candidates, int ele, int *res){
	int count=0;
	int *pd=(int *) malloc(m*sizeof(int));
	getPD(x, m, pd);
	int st=-1;
	dtype *stv = (dtype *) malloc(m*sizeof(dtype));
	int *sti = (int *) malloc(m*sizeof(int));
	for(int e=0;e<ele;e++){
		dtype *s = y+candidates[e];
		st=-1;
		int i;
		for(i=0;i<m;i++){
			dtype value;
			int index;
			int npd;
			while(st!=-1){
				value=stv[st];
				index=sti[st];
				if(value<=s[i]) break;
				st--;
			}
			if(st==-1) pd=0;
			else npd=i-index;
			if(npd!=pd[i]) break;
			stv[++st]=s[i];
			sti[st]=i;
		}
		if(i==m){
			res[count++]=candidates[e];
		}
	}
	return count;
}

#endif
