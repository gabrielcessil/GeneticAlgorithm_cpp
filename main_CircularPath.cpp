#include "genetic.hpp"
#include <cmath>
#include <vector>
#include <chrono>

using namespace std;

double xw = 0, yw = 0;

// Calculate coodinate from the base motors angules
void position(float t1, float t2, float* coordinates){
        double r1 = 5, r2 = 12, r3 = 13.5/2;
        
        double o = 1;
        double e = r1*(sin(t1) - sin(t2)) / ( 2*r3 +r1*cos(t2) -r1*cos(t1) ) ;
        double f = r1*r3*(cos(t2) + cos(t1))/ ( 2*r3 +r1*cos(t2) -r1*cos(t1) );
        double d = 1 + pow(e,2);
        double g = 2*( e*f -e*r1*cos(t1) +e*r3 -r1*sin(t1));
        double h = pow(f,2) -2*f*( r1*cos(t1)-r3 ) -2*r1*r3*cos(t1) + pow(r3,2) + pow(r1,2) - pow(r2,2);

        double y = (-g + o *sqrt( pow(g,2) -4*d*h))/(2*d);
        double x = e*y +f;

        coordinates[0] = x;
        coordinates[1] = y;
}

// User should be aware of the limits for its application: if any gene is out of range the returned value must be OUT_OF_RANGE global;
float coastFunction(array<float, number_of_genes> floatGenes){
        if(	floatGenes[0] > max_gene_value
	     || floatGenes[0] < min_gene_value
	     || floatGenes[1] > max_gene_value
	     || floatGenes[1] < min_gene_value)	return(OUT_OF_RANGE);
        
        
        float coordinate[2];

        position(floatGenes[0], floatGenes[1],  coordinate);

        float x = coordinate[0], y = coordinate[1];

        // The coast is set as the distance from the wanted coordinate to the calculate coordinate
        float coast = sqrt( pow( abs(xw - x), 2 ) + pow( abs(yw - y), 2 ));

        if(isnan(coast)) return(OUT_OF_RANGE);
        return(coast);       
}

// Coordinate calculator for each iteration
float y_fun(float i){
    return(10.5 + 1.5*sin(i));
}
float x_fun(float i){
    return(1.5*cos(i));
}

// CIRCULO 

int main(){
    
    // Save files
    ofstream file_xy, file_angules;
    file_xy.open        ("planningXY.csv", ofstream::trunc | ofstream::in);
    file_angules.open   ("planningAngules.csv", ofstream::trunc | ofstream::in);

    float i = 0;

    while( i < 3.1415*2 ){

        // Define goal point 
        xw = x_fun(i);
        yw = y_fun(i);

        // Auxiliar variables
        float * angules;
        float coordinate[2];

        // Define the genetic algorithm
        GeneticAlgorithm ga(&coastFunction);
        //ga.enableResultDisplay();

        // Run genetic algorithm for the current function coast
        angules = ga.run();

        // Save solution angules
        file_angules << angules[0] << "," << angules[1] << "\n";

        // Save solution coordinate and correct answear
        position(angules[0], angules[1], coordinate);
        file_xy << xw << ", " << yw << ",";
        file_xy << coordinate[0] << "," << coordinate[1] << "\n";

        i += 3.1415*2 / 8;

        cout <<"X: " << coordinate[0] << " Y: " << coordinate[1] << " -> Left motor: " << angules[0] << " Right motor:" << angules[1] << endl;

        delete [] angules;
    
    }
    
    file_xy.close();
    file_angules.close();
    
    return(0);
}
