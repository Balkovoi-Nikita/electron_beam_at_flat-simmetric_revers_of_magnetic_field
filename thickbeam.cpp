#include <stdio.h>
#include <math.h>

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923078164062f
#define Bz(C0, z, x) ((1.0 + C0) / 2.0 + ((1 - C0) / 2.0)*cos(PI*z))
#define Bx(C0, z, x) (PI*((1.0 - C0) / 2.0)*sin(PI*z)*x)

#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

void printProgress(double percentage)
{
	int val = (int)(percentage * 100);
	int lpad = (int)(percentage * PBWIDTH);
	int rpad = PBWIDTH - lpad;
	printf("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
	fflush(stdout);
}
int main()
{
	FILE* params;
	FILE* precise;
	FILE* data;
	FILE* timeweight;
	params  = fopen("params.txt",  "r");
	precise = fopen("precise.txt", "r");
	char filenamedata[60];
	char filenametw[60];

	float C0, nc, w, Xc, Rho;
	fscanf(params, "%f%f%f%f%f", &C0, &nc, &w, &Xc, &Rho);

	sprintf(filenamedata, "Trajectory_C%2.2f_nc%2.2f_w%2.2f_Xc%2.2f_Rho%2.2f.csv", C0, nc, w, Xc, Rho);
	sprintf(filenametw,    "TimeWidth_C%2.2f_nc%2.2f_w%2.2f_Xc%2.2f_Rho%2.2f.csv", C0, nc, w, Xc, Rho);
	printf("%s\n%s\n",filenametw, filenamedata);
	data = fopen(filenamedata, "w");
	timeweight = fopen(filenametw, "w");
	float Rc = sqrt(w) / (PI * 2 * nc);
	Xc  *= Rc;
	Rho *= Rc;

	int eNum, layersNum;
	float Tstep;
	int TrToSkip = 0;
	fscanf(precise, "%d%d%f%d", &eNum, &layersNum, &Tstep, &TrToSkip);
	TrToSkip = eNum / TrToSkip;
	if (TrToSkip == 0)
		TrToSkip = 1;
	int sliceNum = 2 * layersNum + 1;

	printf("C0 = %f\nnc = %f\nRc = %f\nXc = %f\nRho = %f\nTstep = %f\nNumber of electron = %d\n",
		C0, nc, Rc, Xc, Rho, Tstep, eNum*sliceNum);
	fprintf(data, "t, x, y, z, vx, vy, vz, weight\n");
	fprintf(timeweight, "t, weight, xr\n");
	for (int i = -layersNum; i <= layersNum; i++) {
		float xr;
		float weight;
		if (layersNum != 0) {
			xr = (float)i / layersNum;
			weight = sqrt(1 - xr*xr);
			if (i == layersNum || i == -layersNum)
				weight = 0.5f / layersNum;
			//половина засчёт того, что нет слоя сверху
		}
		else {
			xr = 0;
			weight = 1;
		}
		xr *= Rho;
		float t;
		float x;
		float y;
		float z;
		float vx;
		float vy;
		float vz;
		for (int j = 0; j < eNum; j++) {
			 t = 2*PI * j / eNum;
			 x = Rc * cos(t) + Xc + xr;
			 y = Rc * sin(t);
			 z = 0;
			 vx = -Rc * sin(t);
			 vy = Rc * cos(t);
			 vz = 1 / (2*PI * nc);
			while (z < 1 && t < 2 * PI * j / eNum + 2 * PI * nc * 2 && z > -0.01) {
				t += Tstep;

				float k_x1 = vx*Tstep;
				float k_y1 = vy*Tstep;
				float k_z1 = vz*Tstep;
				float k_vx1 = -Bz(C0, z, x)*vy*Tstep;
				float k_vy1 = (Bz(C0, z, x)*vx - Bx(C0, z, x)*vz)*Tstep;
				float k_vz1 = Bx(C0, z, x)*vy*Tstep;

				float k_x2 = (vx + k_vx1/2)*Tstep;
				float k_y2 = (vy + k_vy1/2)*Tstep;
				float k_z2 = (vz + k_vz1/2)*Tstep;
				float k_vx2 = -Bz(C0, z+k_z1/2, x+k_x1/2)*(vy + k_vy1/2)*Tstep;
				float k_vy2 = 
					(Bz(C0, z+k_z1/2, x+k_x1/2)*(vx + k_vx1/2) - Bx(C0, z+k_z1/2, x+k_x1/2)*(vz + k_vz1/2))
					*Tstep;
				float k_vz2 = Bx(C0, z+k_z1/2, x+k_x1/2)*(vy + k_vy1/2)*Tstep;

				float k_x3 = (vx + k_vx2/2)*Tstep;
				float k_y3 = (vy + k_vy2/2)*Tstep;
				float k_z3 = (vz + k_vz2/2)*Tstep;
				float k_vx3 = -Bz(C0, z + k_z2/2, x + k_x2/2)*(vy + k_vy2/2)*Tstep;
				float k_vy3 =
					(Bz(C0, z + k_z2/2, x + k_x2/2)*(vx + k_vx2/2) - Bx(C0, z + k_z2/2, x + k_x2/2)*(vz + k_vz2/2))
					*Tstep;
				float k_vz3 = Bx(C0, z + k_z2/2, x + k_x2/2)*(vy + k_vy2/2)*Tstep;

				float k_x4 = (vx + k_vx3)*Tstep;
				float k_y4 = (vy + k_vy3)*Tstep;
				float k_z4 = (vz + k_vz3)*Tstep;
				float k_vx4 = -Bz(C0, z + k_z3, x + k_x3)*(vy + k_vy3)*Tstep;
				float k_vy4 =
					(Bz(C0, z + k_z3, x + k_x3)*(vx + k_vx3) - Bx(C0, z + k_z3, x + k_x3)*(vz + k_vz3))
					*Tstep;
				float k_vz4 = Bx(C0, z + k_z3, x + k_x3)*(vy + k_vy3)*Tstep;

				x  += k_x1  / 6 + k_x2  / 3 + k_x3  / 3 + k_x4  / 6;
				y  += k_y1  / 6 + k_y2  / 3 + k_y3  / 3 + k_y4  / 6;
				z  += k_z1  / 6 + k_z2  / 3 + k_z3  / 3 + k_z4  / 6;
				vx += k_vx1 / 6 + k_vx2 / 3 + k_vx3 / 3 + k_vx4 / 6;
				vy += k_vy1 / 6 + k_vy2 / 3 + k_vy3 / 3 + k_vy4 / 6;
				vz += k_vz1 / 6 + k_vz2 / 3 + k_vz3 / 3 + k_vz4 / 6;
			if(j%TrToSkip == 0)
				fprintf(data, "%f, %f, %f, %f, %f, %f, %f, %f\n", t, x, y, z, vx, vy, vz, weight);
			}
			fprintf(timeweight, "%f, %f, %f\n", t, weight, xr);
			printProgress((((float)i + layersNum)*eNum + j) / sliceNum / eNum);
		}
	}
	printProgress(1);
	printf("\nDone");
	fclose(precise);
	fclose(data);
	fclose(params);
	fclose(timeweight);
	getchar();
	return 0;
}