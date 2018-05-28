//Uncomment the function to test
#include <stdio.h>
#include "FFbasic.h"

int main(void){
    TestEosPure();
    TestEosMix();
    TestCorrelationOptimization();//takes 20 seconds
    TestEOSoptimization();//takes 30 seconds, or the time you specify in the function
    TestUnifac();

    printf("Press enter to finalize...\n");
    getchar();
    return 0;
}
