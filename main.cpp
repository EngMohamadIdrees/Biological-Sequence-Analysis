#include<cstring>
#include<iostream>
#include<algorithm>
#include<chrono>
#define DIR "D:\\Bioinformatics 3(2)\\Bioligical\\lab 3\\Assignment1\\"
using namespace std;
using namespace std::chrono;
//int Time;
char *Input=new char [10001]; // This array of char is where i will save the whole sequence of size 10000.
char *Input_whole=new char[6000000]; //for whole Genome Sequences
int *suffix_arry=new int [10001]; //This array of Suffixes where i will save all the Suffixes ID.
int *suffix_arry_whole=new int [6000000]; //for whole Genome Sequences
int *Order=new int [10001]; // This array where i will save order of each group of char when i sorted the suffixes
int *Order_whole=new int [6000000]; //for whole Genome Sequences
int *New_Order=new int [10001];//This is a Temporary array when end i will append this array in Order array
int *New_Order_whole=new int [6000000]; // for whole Genome Sequences
char *Input_Test=new char [11];
int *suffix_arry_Test=new int [11];
int *Order_Test=new int [11];
int *New_Order_Test=new int [11];


int Curent_len;
char * Read(char *Input)
{
    char *buf=new char[10001];   // This buffer is used to Temporary store the sequence line by line
    FILE* file=fopen(DIR "genome.fastq", "r");
    int Begin=0;  // This variable is used to know which line i will be in
    int Start=0;  // This variable is used to know where i am reach in the array of char
    int Stop=0;   // This variable is used to know when i will stop the read form the file
    while(!feof(file)&&Stop!=1) // this while loop will loop in the file input and it will not stop until the Input size of sequence reach 10000
    {

        if(Begin==0) //here i will check if i will be in the first line which have the id of the gene sequences
        {
            fscanf(file, "%[^\n\r] ", buf);
            buf[0]=0;
            Begin++;
            continue;   // i will discard this line be read this sequence in buffer then i will clear the buffer and continue to read other lines
        }
        else
        {
            fscanf(file, "%[^\n\r] ", buf);
            for(int i=0; i<strlen(buf); i++)
            {
                if(strlen(Input)==10000) // here i check if the size of Input reach 10000 that means i have now 10000 char of sequence i will let 1 char for the sentinel $
                {
                    Input[Start]='$';   // I add the $ depending on  that the Input array is have sequence of size 10001
                    Stop=1; //here i will stop the while loop
                    break;
                }
                Input[Start]=buf[i]; // if the size of Input did not reach the 10000 that mean i will continue read of char of the sequence gene in the buf
                Start++;
            }
        }
    }
    Input[Start]='$';      // This is alternatives above,if we want to add the ($).
    fclose(file);
    fflush(NULL);
    delete[]buf; // here i will destroy the buffer
    return Input;
}
char * Read_whole(char *Input)
{
    char *buf_whole=new char[6000000];   // This buffer is used to Temporary store the sequence line by line
    FILE* file=fopen(DIR "genome.fastq", "r");
    int Begin=0;  // This variable is used to know which line i will be in
    int Start=0;  // This variable is used to know where i am reach in the array of char
    int Stop=0;   // This variable is used to know when i will stop the read form the file
    while(!feof(file)&&Stop!=1) // this while loop will loop in the file input and it will not stop until the Input size of sequence reach the end of sequences
    {

        if(Begin==0) //here i will check if i will be in the first line which have the id of the gene sequences
        {
            fscanf(file, "%[^\n\r] ", buf_whole);
            buf_whole[0]=0;
            Begin++;
            continue;   // i will discard this line be read this sequence in buffer then i will clear the buffer and continue to read other lines
        }
        else
        {
            fscanf(file, "%[^\n\r] ", buf_whole);
            for(int i=0; i<strlen(buf_whole); i++)
            {
                if(buf_whole[i]=='>') //here this is my base when i reach the > i will stop read form file
                {
                    Stop=1;
                    break;
                }
                Input[Start]=buf_whole[i];
                Start++;
            }
        }
    }
    Input[Start]='$';      // This is alternatives above,if we want to add the ($).
    fclose(file);
    fflush(NULL);
    delete[]buf_whole; // here i will destroy the buffer
    return Input;
}
char * Read_Test(char *Input,int Line)
{
    char *buf_Test=new char[90];   // This buffer is used to Temporary store the sequence line by line
    FILE* file=fopen(DIR "genome.fastq", "r");
    int Begin=0;  // This variable is used to know which line i will be in
    int Start=0;  // This variable is used to know where i am reach in the array of char
    int Stop=0;   // This variable is used to know when i will stop the read form the file
    while(!feof(file)&&Stop!=1) // this while loop will loop in the file input and it will not stop until the Input size of sequence reach the end of sequences
    {

        while(Begin<Line) //here i will check if i will be in the first line which have the id of the gene sequences
        {
            fscanf(file, "%[^\n\r] ", buf_Test);
            buf_Test[0]=0;
            Begin++;
            continue;   // i will discard this line be read this sequence in buffer then i will clear the buffer and continue to read other lines
        }
        fscanf(file, "%[^\n\r] ", buf_Test);
        //cout<<strlen(buf_Test)<<endl;
        for(int i=0; i<strlen(buf_Test); i++)
        {
            if(i==10) //here this is my base when i reach the > i will stop read form file
            {
                Stop=1;
                break;
            }
            Input[Start]=buf_Test[i];
            Start++;
        }
    }
    Input[Start]='$';      // This is alternatives above,if we want to add the ($).
    fclose(file);
    fflush(NULL);
    delete[]buf_Test; // here i will destroy the buffer
    return Input;
}

void print(char *C_Input) // in this function i can print all the sequences i had read it
{
    for(int i=0; i<strlen(C_Input); i++)
    {
        cout<<C_Input[i];
    }
    cout<<endl;
}
int *build_suffix(char* Current_Input,int *Current_suffix_arry) //in this function i will expand all the suffix of the Input sequences
{
    int Size=strlen(Current_Input);
    for(int i=0; i<Size; i++)
    {
        Current_suffix_arry[i]=i;
    }
    return Current_suffix_arry;   // here this array will have all the  start index of each suffixes
}
void Print_Suffixes(int * Current_suffix_arry,char *Current_Input)  // in this function i will print all the Suffixes with each ID
{
    for(int i=0; i<strlen(Current_Input); i++)
    {
        cout<<"Suffix ID ="<<Current_suffix_arry[i]<<endl;
        for(int j=Current_suffix_arry[i]; j<strlen(Current_Input); j++)
        {
            cout<<Current_Input[j];
        }
        cout<<endl;
        cout<<endl;
    }
}
bool Compare(int a,int b,char *Current_Input)
{
    int Size=strlen(Current_Input);
    while(a<Size&&Current_Input[a]==Current_Input[b]) // Compare the two Suffixes until find the smallest one of them  or when the length of one of them  has ended
    {
        if (b>=Size) break;
        ++a;
        ++b;

    }
    return Current_Input[a]<Current_Input[b];
}
void merge(int* arr, int l, int m, int r,char *Current_Input)// in this function i will merge and sort the suffixes depending on Alphabet
{
    int i, j, k;
    int n1 = m - l + 1;// First subsequence array will be form l to  m  --> arr[l..m]
    int n2 = r - m;// Second subsequence array  will be form m+1 to r ---> arr[m+1..r]
    int *L=new int[n1], *R=new int[n2]; // i will create two arrays to merge in

    for (i = 0; i < n1; i++)
        L[i] = arr[l + i];
    for (j = 0; j < n2; j++)
        R[j] = arr[m + 1 + j];

    i = j = 0;
    k = l;

    while (i < n1 && j < n2)
    {
        if (Compare(L[i],R[j],Current_Input)) //here i merged the Suffixes depending on Alphabet
        {
            arr[k] = L[i];
            i++;
        }
        else
        {
            arr[k] = R[j];
            j++;
        }
        k++;
    }

    while (i < n1) //this while is used to if the on of the two Suffixes is Long than one
    {
        arr[k] = L[i];
        i++;
        k++;
    }

    while (j < n2)  //this while is used to if the on of the two Suffixes is Long than one
    {
        arr[k] = R[j];
        j++;
        k++;
    }
}
void mergeSort(int* arr, int l, int r,char *Current_Input) //here this function is mergeSort Recursion will be  run recursively until have the sorted Suffixes array
{
    if (l < r)
    {
        int m = l + (r - l) / 2;

        mergeSort(arr, l, m,Current_Input);
        mergeSort(arr, m + 1, r,Current_Input);

        merge(arr, l, m, r,Current_Input);
    }
}
bool Compare2(int a,int b) //here i will compare the two group of char have two ID depending on the ASCI of each char
{
    return Order[a]<Order[b]|| Order[a]==Order[b]&&Order[a+Curent_len]<Order[b+Curent_len];
}
void predfix_doubling(int *suffix_arry)
{
    int n=strlen(Input);
    for(int i=0; i<n; i++)
    {
        Order[i]=Input[i]; // i will save all the ASCI of each char of the Input Sequences
    }
    memset(New_Order,0,n); // i will create temporary array with size n to save the orders in each step
    for (int Len=1;; Len*=2)
    {
        Curent_len=Len;
        sort(suffix_arry,suffix_arry+n,Compare2); //First i will sort depend on the first char then i will sort depend on the two char and then will keep in 2^n iteration
        for(int i=1; i<strlen(Input); i++)
        {
            New_Order[i]=New_Order[i-1]+Compare2(suffix_arry[i-1],suffix_arry[i]); //the Order of the current is depending on the result of compare, will be == the last or greater than bay one
        }
        for(int i=0; i<strlen(Input); i++)
        {
            Order[suffix_arry[i]]=New_Order[i]; // in the end of iteration i will save the order of each suffixe
        }
        if(New_Order[n-1]==n-1)
        {
            break;
        }

    }
}
bool Compare2_whole(int a,int b) //here i will compare the two group of char have two ID depending on the ASCI of each char
{
    return Order_whole[a]<Order_whole[b] || Order_whole[a]==Order_whole[b]&&Order_whole[a+Curent_len]<Order_whole[b+Curent_len];
}
void predfix_doubling_whole(int *suffix_arry_whole)
{
    int n=strlen(Input_whole);
    int z=0;
    for(int i=0; i<n; i++)
    {
        Order_whole[i]=Input_whole[i]; // i will save all the ASCI of each char of the Input Sequences
    }
    memset(New_Order_whole,0,n); // i will create temporary array with size n to save the orders in each step
    for (int Len=1;; Len*=2)
    {
        Curent_len=Len;
        sort(suffix_arry_whole,suffix_arry_whole+n,Compare2_whole); //First i will sort depend on the first char then i will sort depend on the two char and then will keep in 2^n iteration
        for(int i=1; i<n; ++i)
        {
            New_Order_whole[i]=New_Order_whole[i-1]+Compare2_whole(suffix_arry_whole[i-1],suffix_arry_whole[i]); //the Order of the current is depending on the result of compare, will be == the last or greater than bay one
        }
        for(int i=0; i<n; ++i)
        {
            Order_whole[suffix_arry_whole[i]]=New_Order_whole[i]; // in the end of iteration i will save the order of each suffixe
        }
        if(New_Order_whole[n-1]==n-1)
        {
            break;
        }

    }
}
void Print__Suffixes_One_In_file(int *Suffixes,long long Time)
{
    int Length=strlen(Input);
    FILE *in_file  = fopen("Suffixes_1.txt", "w"); // Open new file
    fprintf(in_file,"Sorted Suffixes_1 with Complexity = O(n^2log(n))\n"); // print first line in the file
    fprintf(in_file,"Time taken in Sort Function by Merge in Mircroseconds = %d \n",Time); // print second line in the file
    for(int i=0; i<Length; i++)
    {
        fprintf(in_file,"Suffix ID: %d\n",Suffixes[i]); // print the suffixes ID's
    }

    fclose(in_file);
}
void Print__Suffixes_Two_In_file(int *Suffixes,long long Time)
{
    int Length=strlen(Input);
    FILE *in_file  = fopen("Suffixes_2.txt", "w"); // Open new file
    fprintf(in_file,"Sorted Suffixes_2 with Complexity = O(n log^2 n)\n"); // print first line in the file
    fprintf(in_file,"Time taken in Sort Function by Predifx in Mircroseconds = %d \n",Time); // print second line in the file
    for(int i=0; i<Length; i++)
    {
        fprintf(in_file,"Suffix ID: %d\n",Suffixes[i]); // print the suffixes ID's
    }
    fclose(in_file);
}
void Print__Suffixes_Whole_In_file(int *Suffixes,long long Time)
{
    int Length=strlen(Input_whole);
    FILE *in_file  = fopen("Suffixes_3.txt", "w"); // Open new file
    fprintf(in_file,"Sorted Suffixes_3 with Complexity = O(n log^2 n)\n"); // print first line in the file
    fprintf(in_file,"Time taken in Sort Function by predfix in seconds = %d \n",Time); // print second line in the file
    for(int i=0; i<Length; i++)
    {
        fprintf(in_file,"Suffix ID: %d\n",Suffixes[i]); // print the suffixes ID's
    }
    fclose(in_file);
}
bool Compare2_Test(int a,int b) //here i will compare the two group of char have two ID depending on the ASCI of each char
{
    return Order_Test[a]<Order_Test[b]|| Order_Test[a]==Order_Test[b]&&Order_Test[a+Curent_len]<Order_Test[b+Curent_len];
}
void predfix_doubling_test(int *suffix_arry_Test)
{
    int n=strlen(Input_Test);
    for(int i=0; i<n; i++)
    {
        Order_Test[i]=Input_Test[i]; // i will save all the ASCI of each char of the Input Sequences
    }
    memset(New_Order_Test,0,n); // i will create temporary array with size n to save the orders in each step
    for (int Len=1;; Len*=2)
    {
        Curent_len=Len;
        sort(suffix_arry_Test,suffix_arry_Test+n,Compare2_Test); //First i will sort depend on the first char then i will sort depend on the two char and then will keep in 2^n iteration
        for(int i=1; i<n; i++)
        {
            New_Order_Test[i]=New_Order_Test[i-1]+Compare2_Test(suffix_arry_Test[i-1],suffix_arry_Test[i]); //the Order of the current is depending on the result of compare, will be == the last or greater than bay one
        }
        for(int i=0; i<n; i++)
        {
            Order_Test[suffix_arry_Test[i]]=New_Order_Test[i]; // in the end of iteration i will save the order of each suffixe
        }
        if(New_Order_Test[n-1]==n-1)
        {
            break;
        }

    }
}

void Delete(char * Sequence, int * Suffixes) // Finally when i done i call this function to destroy both the Input and Suffixes arrays
{
    delete [] Sequence;
    delete [] Suffixes;
}
int main()
{
    // **********here i will do sort in the 10,000 sequences suffixes once using merge sort and once suing prifixs doubling********

    //Input=Read(Input); // with this function i will read 10000 of Gene Sequences.

    //print(Input);     //here i can print all the Sequences i had read it.

    //cout<< strlen(Input)<<endl;  //here i can print the length of my sequences.

    //suffix_arry=build_suffix(Input,suffix_arry);     //here i build the suffix array.

    //cout<<"Before sort ="<<endl;

    //Print_Suffixes(suffix_arry,Input);    //this print will print the suffixes before sort it.

    //auto start = std::chrono::system_clock::now(); //here where i will star calculate the time before using merge sort

    //mergeSort(suffix_arry,0,strlen(Input)-1,Input);  //here i will sort the suffixes depending on the Alphabet in time O(n^2 log n).

    //auto end = std::chrono::high_resolution_clock::now()-start; //here i will stop the time after i sorted the array

    //long long Time=std::chrono::duration_cast<std::chrono::microseconds>(end).count(); // i calculate the time by subtract the end from start

    //cout<<"After Sort ="<<endl;

    //Print_Suffixes(suffix_arry,Input);      //this print will print the suffixes after i sort it.

    //Print__Suffixes_One_In_file(suffix_arry,Time); // this print will print in file suffixes_1 all sorted suffixes with each id.

    //auto start = std::chrono::high_resolution_clock::now();

    //predfix_doubling(suffix_arry);    //here i will sort the suffixes in time O(n log^2 n).

    //auto end = std::chrono::high_resolution_clock::now()-start;

    //long long Time=std::chrono::duration_cast<std::chrono::microseconds>(end).count();


    //Print_Suffixes(suffix_arry,Input);        //this print will print the suffixes after i sort it.

    //Print__Suffixes_Two_In_file(suffix_arry,Time); // this print will print in file suffixes_2 all sorted suffixes with each id.

    //Delete(Input,suffix_arry);      // Finally i delete and destroy the both of Sequences and Suffixes

    //--------------------------------------------------------------------------------------------------------------------------------

     // ******************here i will do sort in the whole sequences suffixes using prifix doubling************************

    //Input_whole=Read_whole(Input_whole);

    //print(Input_whole);

    //suffix_arry_whole=build_suffix(Input_whole,suffix_arry_whole);

    //cout<<"Before sort ="<<endl;

    //Print_Suffixes(suffix_arry_whole,Input_whole);

    //auto start = std::chrono::high_resolution_clock::now();

    //predfix_doubling_whole(suffix_arry_whole);    //here i will sort the suffixes in time O(n log^2 n).

    //Print_Suffixes(suffix_arry,Input);        //this print will print the suffixes after i sort it.

    //auto end = std::chrono::high_resolution_clock::now()-start;

    //long long Time=std::chrono::duration_cast<std::chrono::microseconds>(end).count();

    //Print__Suffixes_Whole_In_file(suffix_arry_whole,Time); // this print will print in file suffixes_3 all sorted suffixes with each id.

    //Delete(Input,suffix_arry);      // Finally i delete and destroy the both of Sequences and Suffixes.

    //----------------------------------------------------------------------------------------------------------------------------------

    // Test Case :
    /*
    for(int i=1;i<=30;i++)
    {
        Input_Test=Read_Test(Input_Test,i);

        cout<<"Test Case "<<i<<" : "<<endl;


        cout<<endl;

        cout<<"The sequence is : ";

        cout<<Input_Test<<endl;

        cout<<endl;

        suffix_arry_Test=build_suffix(Input_Test,suffix_arry_Test);

        cout<<"Before Sort the Suffixes :"<<endl;

        cout<<endl;

        Print_Suffixes(suffix_arry_Test,Input_Test);

        mergeSort(suffix_arry_Test,0,strlen(Input_Test)-1,Input_Test);

        cout<<"After sort the suffixes by using Merge sort :"<<endl;

        cout<<endl;

        Print_Suffixes(suffix_arry_Test,Input_Test);

        predfix_doubling_test(suffix_arry_Test);

        cout<<"After sort the suffixes by using predfix sort :"<<endl;

        cout<<endl;

        Print_Suffixes(suffix_arry_Test,Input_Test);

        Delete(Input_Test,suffix_arry_Test);




    }
    */
    return 0;
}
