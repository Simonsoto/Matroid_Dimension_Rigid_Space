'''Usefull information of matroid functions already implemented
.circuits() Return the collection of circuits.
.cocircuits() Return the collection of cocircuits.
.nonbases() Set of nonbases of a matroid.
'''


def finding_circuit(non_basis, circuits):
    '''
    This function find the circuit that is inside the non-basis of rank r-1
    Input:
    
    circuits: List of Sets that represent circuits of my matroid. 
    non_basis: The Set representing the non_basis of the matroid that we want to study.
    
    return:
    the circuit that is inside the non-basis.
    '''
    for circuit in circuits:
        if Set(circuit) in non_basis.subsets():
            return Set(circuit)
        

def finding_cocircuit(non_basis, cocircuits, ground_set):
    '''
    This function find the cocircuit that is inside the non-basis of rank r-1
    Input:
    
    cocircuits: List of sets that represent circuits of my matroid.
    non_basis: The non_basis of the matroid that we want to study.
    ground_set: Ground set of the matroid.
    
    return:
    the cocircuit that is inside the non-basis.
    '''
    complement=ground_set.symmetric_difference(non_basis)
    for cocircuit in cocircuits:
        if Set(cocircuit) in complement.subsets():
            return Set(cocircuit)
#Now we can calculate the symbols associated to the corresponding non_basis and the circuit and cocircuit associated.
def cal_symbols_non_basis(non_basis, circuit, cocircuit, symbols):
    '''
    Function to calculate the collection of symbols of a given non_basis.
    
    Input:
    non_basis: Non_basis of the matroid, represented as a set from which we are calculating the symbols
    circuit: Set representing the circuit contained in the non_basis
    cocircuit: Set representing the cocircuit contained in the complement of the non_basis.
    symbols: Collection of symbols that satisfies the required condition.
    '''
    #Particular element circuit
    a=Set(circuit[0])
    rest_circuit=circuit.difference(a)
    
    #particular element cocircuit
    d=Set(cocircuit[0])
    rest_cocircuit=cocircuit.difference(d)
        
    for b in rest_circuit:
        b=Set(b)
        S=non_basis.difference(a.union(b))        
        
        for c in rest_cocircuit:
            c=Set(c)
#             print([S,a.union(b),d.union(c)])
            symbols.append([S,a.union(b),d.union(c)])
    return symbols

#Find the set of desired non-bases to be used
def non_bases_ideal_rank(matroid,rank,non_bases):
    '''Function to find the collection of non_basis with rank r-1 where r is the rank of the matroid.
    Input:
    matroid: the matroid that we are working with
    rank: Integer that represents the rank ofthe matroid
    non_bases: Collection of sets that represent the non-bases of the matroid.
    
    return:
    new collection of non-basis that all of them have rank r-1
    '''
    special_non_bases=[]
    for non_basis in non_bases:
        if matroid.rank(non_basis)==rank-1:
            special_non_bases.append(non_basis)
    return special_non_bases
    
    
#Vector space for each set of symbols! 
#The idea for this is that each basis is a basis of the vector space RR^B
#In this vector space

#First create the set of vectors associated to a collection of symbols.
def equation_vectors(bases, symbols, dim_B):
    '''
    Function that calculate the vectors associated with a set of symbols of a non-basis..
    Input:
    bases: Collection of sets that represent the bases of the matroid. Give the order for the vector space.
    symbols_non_basis: Symbols associated.
    dim_B: Dimension of the ambient vector space. Number of bases.
    Output:
    vectors: Collection of vectors that represent the equations given by the symbols.
    '''
    
    vectors=[]
    for symbol in symbols:
        vector=[0]*dim_B
        #extracting the information out of the symbol
        S=symbol[0]
        a=Set(symbol[1][0])
        b=Set(symbol[1][1])
        c=Set(symbol[2][0])
        d=Set(symbol[2][1])
        #Creating the four basis associated to the symbol.
        Sac=S.union(a).union(c)
        Sad=S.union(a).union(d)
        Sbc=S.union(b).union(c)
        Sbd=S.union(b).union(d)
        # Finding the indices of the different basis in the list of bases
        
        #Changing the position in vector positive numbers
        iSac=bases.index(Sac)
        vector[iSac]=1
        iSbd=bases.index(Sbd)
        vector[iSbd]=1
        
        #Changing the position in vector negative numbers
        iSad=bases.index(Sad)
        vector[iSad]=-1
        iSbc=bases.index(Sbc)
        vector[iSbc]=-1
        
        #Add the current vector to the list.
        vectors.append(vector)
        
    return vectors

'''
Start grouping the functions already created to calculate the final matrix.
'''
def equations_vector_space(matroid,bases,dim_B,ground_set):
    '''
    Function that calculate all the vectors in the vector space associated to the matroid.
    Input:
    matroid: Matroid. The matroid that the vectors are going to be calculated.
    bases: List of Sets. Bases of the matroid.
    dim_B: Integer. Represents the amoun of bases of the matroid.
    Output:
    vectors: Collection of vectors that represent the equations of the hyperplanes defined by the 
    special U(2,4) minors of the matroid
    
    '''
    # Basic Matroid information

    # rank of the matroid
    rank=matroid.rank()
    #non bases of my matroid
    non_bases=list(matroid.nonbases())
    #Special set of non bases that have rank rank(M)-1.
    important_non_bases=non_bases_ideal_rank(matroid, rank, non_bases)
    
    #Collection of circuits of the matroid
    circuits=list(matroid.circuits())
    #Collection of cocircuits of the matroid
    cocircuits=list(matroid.cocircuits())
    
    
    symbols=[]
    #Calculate the vectors associated to the collection of special non bases
    for non_basis in important_non_bases:
        non_basis=Set(non_basis)
        #finding circuit contained in non basis
        cir=finding_circuit(non_basis, circuits)
        #finding cocircuit contained in the complement of the non basis
        coc=finding_cocircuit(non_basis, cocircuits, ground_set)
        symbols=cal_symbols_non_basis(non_basis, cir, coc, symbols)
    
    #Check the uniqueness of the elements in the collection of symbols
    
    
    #Initialize collection of vectors for the vector space.
    vectors=equation_vectors(bases,symbols,dim_B )
    
    return vectors

def dimension_D_M(matroid):
    
    
    #bases of the matroid
    bases=list(matroid.bases())
    #Dimension of vector space the vectors
    dim_B=len(bases)
    #Size of the ground set
    #ground set of the matroid first
    ground_set=Set(matroid.groundset())
    #Cardinality of the set.
    size=ground_set.cardinality()
    
    #Finding the equation of the vectors associated to the matrix.
    vectors=equations_vector_space(matroid, bases, dim_B,ground_set)
    
    #uniqueness of the vectors
    unique=[]
    for vector in vectors:
        if vector not in unique:
            unique.append(vector)
    
    #Creating the matrix.
    matrix=Matrix(vectors)
    #Transposing the matrix to be able to calculate the kernel of it
    #In sage ker(matrix) is the left kernel of a amtrix.(x*A=0)
    matrix=matrix.transpose()
    
    #d(matroid) as defined in the talk.
    dimension = dim(kernel(matrix))-size
    return dimension


'''
Keys of the dictionary of matroids up to 9 elements:
Comments on the file matroid_9.sobj
It is a dictionary with the following format
key: a pair(a,b) with a+b<=9 and a is the rank of the matroids, a+b is the number of elements on the matroid
Value: List of matroids that satisfy the characteristics of the key.


'''




''' 
Function that given a matroid on up to nine elements, change the ground set to work with strings instead
of integers as they may cause problems with the other code. 

Input: Matroid. Ground set can be any type of set.
return: A matroid with a suitable ground set for the other code. 
'''
def transformate_matroid_to_letter(matroid):
    #Letters to have as a groundset
    alphabet='abcdefghi'
    #Ground set of the matroid
    ground_set=Set(matroid.groundset())
    #Cardinality of the set
    l=ground_set.cardinality()
    ground_set_copy=Set()
    
    #Make a ground set based on letters
    for let in alphabet[0:l]:
        ground_set_copy=ground_set_copy.union(Set(let))
    
    #create a dictionary with the indicator of the transformation
    dic={}
    for i in range(0,l):
        dic[ground_set[i]]=ground_set_copy[i]
        
    #Change the collection of basis of the matroid for the one in the dictionary
    bases=list(matroid.bases())
    bases_copy=[]
    for basis in bases:
        copy=Set(basis)
        to_add=Set()
        for letter in copy:
            to_add=to_add.union(Set(dic[letter]))
        bases_copy.append(to_add)
        
    #Define a matroid with underlying set as string.
    matroid_copy=Matroid(ground_set_copy,bases_copy)
    
    return matroid_copy
    


    

'''
Function that take the list of matroids on up to 9 elements and calculate the dimension fo D_M
as presented in previously.
Input: Dictionary with the matroids up to nine elements

Output: Dictionary with the matroids and their respective value.
'''
def matroids_up_to_nine(matroid_dict):
    #Making the copy that we are going to return
    matroid_copy=matroid_dict.copy()
    print('making the copy')
    #The keys from the dictionary
    keys=matroid_dict.keys()
        
    
    #working with the list of matroids that is in every key.
    print('working with each of the keys that have at least 4 elements and at least rank 2')
    for key in keys:
        #For each key of the dictionary, we extract the required information.
        r=key[0]
        
        n=key[0]+key[1]
        
        #Algorithm make sense for matroid of rank at least 2 and with at least 4 elements.
        if (r>=2 and n>=4):
            
            #list of matroids of rank r on n elements
            matroids=matroid_dict[key]
            #Work with each matroid.
            #list to put in the dictionary key that we are looking
            new_info=[]
            for element in range(len(matroids)):
                #Make sure that we can work with it by changing the ground set to the one we created the programs for
                matroid= transformate_matroid_to_letter(matroids[element])
                new_info.append([matroids[element],dimension_D_M(matroid)])
                
                
                
            matroid_copy[key]=new_info
        print('New information added')
    return matroid_copy

matroids= load('matroid_9.sobj')
matroids_up_to_nine(matroids)