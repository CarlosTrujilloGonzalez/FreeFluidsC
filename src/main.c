//Uncomment the function to test
int main(void){
    TestEosPure();
    TestEosMix();
    CorrelationOptimization();//takes 20 seconds
    EOSoptimization();//takes 30 seconds, or the time you specify in the function
    printf("Press enter to finalize...\n");
    getchar();
    return 0;
}
