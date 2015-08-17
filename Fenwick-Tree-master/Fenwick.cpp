//Demonstrating the operations of Binary Index Tree
#include <iostream>
using namespace std;

int getTotal(int B_I_Tree[], int n, int index);

void updateB_I_T(int *BITree, int n, int index, int val);

int *constructB_I_Tree(int arr[], int n);

int main()
{
    int arr[] = {2, 1, 1, 3, 2, 3, 4, 5, 6, 7, 8, 9};
    int n = sizeof(arr)/sizeof(arr[0]);
    int *B_I_Tree = constructB_I_Tree(arr, n);
    cout << "The sum of elements in arr[0..5] is "<< getTotal(B_I_Tree, n, 5);
    
    // Testing the update operation
    arr[3] += 6;
    updateB_I_T(B_I_Tree, n, 3, 6); //Updating BIT for above change in arr[]
    
    cout << "\n Sum of the elements in arr[0..5] after update is "<< getTotal(B_I_Tree, n, 5);
    
    return 0;
}
int *constructB_I_Tree(int arr[], int n)
{
    // Initialize B_I_Tree[] as 0
    int *B_I_Tree = new int[n+1];
    for (int i=1; i<=n; i++)
        B_I_Tree[i] = 0;
    
    // Store the actual values in B_I_Tree[]
    for (int i=0; i<n; i++)
        updateB_I_T(B_I_Tree, n, i, arr[i]);
    
    return B_I_Tree;
}
int getTotal(int B_I_Tree[], int n, int index)
{
    int total = 0;
    
    // index in B_I_Tree[] is 1 more than the index in arr[]
    index = index + 1;
    
    while (index>0)
    {
        total += B_I_Tree[index];
        // Move index to parent node
        index -= index & (-index);
    }
    return total;
}
void updateB_I_T(int *B_I_Tree, int n, int index, int val)
{
    // index in B_I_Tree[] is 1 more than the index in arr[]
    index = index + 1;
    
    // Traverse all ancestors and add 'val'
    while (index <= n)
    {
        // Add 'val' to current node of B_I_Tree
        B_I_Tree[index] += val;
        
        // Update index to that of parent
        index += index & (-index);
    }
}
