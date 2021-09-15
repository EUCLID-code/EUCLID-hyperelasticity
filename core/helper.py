from .__importList__ import *

def initCUDA(cudaID):
    """
    
    Setup GPU-CUDA acceleration.
    
    _Input Arguments_
    
    - `cudaID` - -1: CPU, 0: CUDA0, 1: CUDA1
    
    _Output Arguments_
    
    - _none_
    
    ---
    
    """
    print('\n-----------------------------------------------------')
    if(cudaID >= 0):
        device = torch.device('cuda:'+str(cudaID)) #cuda
        torch.cuda.set_device(cudaID) #cuda
        torch.set_default_tensor_type(torch.cuda.DoubleTensor) #cuda
    else:
        device = torch.device('cpu')
        torch.set_default_tensor_type(torch.DoubleTensor)
    print("Setting device to: ",cudaID)
    x=torch.rand(1,1)
    print('Test: ',x.device)
    print('-----------------------------------------------------\n')



