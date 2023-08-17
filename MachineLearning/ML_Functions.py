"""
Functions for the Mathine Learning module

Created on 2023/3/25

@author: Pengfei Jin
"""

import numpy as np;  
import joblib
import time
import concurrent.futures
import math
import torch
from collections import Counter
from torch.utils.data import Dataset, DataLoader, SubsetRandomSampler
import torch.nn as nn



def LinearReg(X_all,Y_all,ind_used,if_bias=True,lambda_ls=1,RoundHalf=False):
    n_total = X_all.shape[0]
    if(if_bias==False):
        X_t = X_all[:,ind_used]
    else:
        X_t = X_all[:,ind_used]
        X_t = torch.cat([X_all,torch.ones([n_total,1])],dim=1)
    n_feature = X_t.shape[1]
    I_hat = torch.eye(n_feature)
    if(if_bias==True):
        I_hat[n_feature-1,n_feature-1] = 0
    m = (torch.matmul(X_t.t(),X_t)+lambda_ls*I_hat)
    m_inv = torch.inverse(m)
    alpha_t = torch.matmul(m_inv,torch.matmul(X_t.t(),Y_all))
    print('Coefficient = ' + ', '.join(['{:.2f}'.format(i) for i in torch.round(alpha_t * 100) / 100]))
    pre = torch.matmul(X_t,alpha_t)
    if RoundHalf==False:
        pre_r = torch.round(pre)
        c_n = (pre_r==Y_all).sum()
    else: 
        pre_r = torch.round(pre*2.0)
        Y_2 = Y_all*2
        c_n = (pre_r==Y_2).sum()
    error = pre-Y_all
    error_all = error.abs().sum()
    print('total: %d'%(n_total))
    print('acc: %.2f'%(c_n.float()/n_total))
    print('average error: %.2f'%(error_all/n_total))
    return



def Dataset_spit(X,Y,ratio=0.8):
    indices = torch.randperm(X.shape[0])
    split_point = int(ratio * X.shape[0])
    X_tr = X[indices[:split_point]]
    X_te = X[indices[split_point:]]
    Y_tr = Y[indices[:split_point]]
    Y_te = Y[indices[split_point:]]
    return X_tr,X_te,Y_tr,Y_te


class LinearF(nn.Module):
    def __init__(self, c_in, c_out):
        super(LinearF, self).__init__()
        self.linear1 = nn.Linear(c_in,c_out)
    def forward(self, x):
        x = self.linear1(x)
        return x


class OneLayerReLU(nn.Module):
    def __init__(self, c_in, c_hidden, c_out):
        super(OneLayerReLU, self).__init__()
        self.linear1 = nn.Linear(c_in,c_hidden)
        self.linear2 = nn.Linear(c_hidden,c_out)
    def forward(self, x):
        x = torch.nn.functional.relu(self.linear1(x))
        x = self.linear2(x)
        return x

class TwoLayerReLU(nn.Module):
    def __init__(self, c_in, c_hidden1, c_hidden2, c_out):
        super(TwoLayerReLU, self).__init__()
        self.linear1 = nn.Linear(c_in,c_hidden1)
        self.linear2 = nn.Linear(c_hidden1,c_hidden2)
        self.linear3 = nn.Linear(c_hidden2,c_out)
    def forward(self, x):
        x = torch.nn.functional.relu(self.linear1(x))
        x = torch.nn.functional.relu(self.linear2(x))
        x = self.linear3(x)
        return x
    
class ThreeLayerReLU(nn.Module):
    def __init__(self, c_in, c_hidden1, c_hidden2, c_hidden3, c_out):
        super(ThreeLayerReLU, self).__init__()
        self.linear1 = nn.Linear(c_in,c_hidden1)
        self.linear2 = nn.Linear(c_hidden1,c_hidden2)
        self.linear3 = nn.Linear(c_hidden2,c_hidden3)
        self.linear4 = nn.Linear(c_hidden3,c_out)
    def forward(self, x):
        x = torch.nn.functional.relu(self.linear1(x))
        x = torch.nn.functional.relu(self.linear2(x))
        x = torch.nn.functional.relu(self.linear3(x))
        x = self.linear4(x)
        return x

class FourLayerReLU(nn.Module):
    def __init__(self, c_in, c_hidden1, c_hidden2, c_hidden3, c_hidden4, c_out):
        super(FourLayerReLU, self).__init__()
        self.linear1 = nn.Linear(c_in,c_hidden1)
        self.linear2 = nn.Linear(c_hidden1,c_hidden2)
        self.linear3 = nn.Linear(c_hidden2,c_hidden3)
        self.linear4 = nn.Linear(c_hidden3,c_hidden4)
        self.linear5 = nn.Linear(c_hidden4,c_out)
    def forward(self, x):
        x = torch.nn.functional.relu(self.linear1(x))
        x = torch.nn.functional.relu(self.linear2(x))
        x = torch.nn.functional.relu(self.linear3(x))
        x = torch.nn.functional.relu(self.linear4(x))
        x = self.linear5(x)
        return x
    
class FiveLayerReLU(nn.Module):
    def __init__(self, c_in, c_hidden1, c_hidden2, c_hidden3, c_hidden4, c_hidden5, c_out):
        super(FiveLayerReLU, self).__init__()
        self.linear1 = nn.Linear(c_in,c_hidden1)
        self.linear2 = nn.Linear(c_hidden1,c_hidden2)
        self.linear3 = nn.Linear(c_hidden2,c_hidden3)
        self.linear4 = nn.Linear(c_hidden3,c_hidden4)
        self.linear5 = nn.Linear(c_hidden4,c_hidden5)
        self.linear6 = nn.Linear(c_hidden5,c_out)
    def forward(self, x):
        x = torch.nn.functional.relu(self.linear1(x))
        x = torch.nn.functional.relu(self.linear2(x))
        x = torch.nn.functional.relu(self.linear3(x))
        x = torch.nn.functional.relu(self.linear4(x))
        x = torch.nn.functional.relu(self.linear5(x))
        x = self.linear6(x)
        return x    
    
    
def NN_train(X,Y,n_layer=1,c_mid=10,learning_rate=0.001,weight_decay=0.01,N_epoch=1,batch_size=32,opt_adam=True,print_mid_loss=False):
    N = X.shape[0]
    c_in = X.shape[1]
    if n_layer==1:
        net = OneLayerReLU(c_in,c_mid,1)
    elif n_layer==2:
        net = TwoLayerReLU(c_in,c_mid,c_mid,1)
    elif n_layer==3:
        net = ThreeLayerReLU(c_in,c_mid,c_mid,c_mid,1)
    elif n_layer==4:
        net = FourLayerReLU(c_in,c_mid,c_mid,c_mid,c_mid,1)
    elif n_layer==5:
        net = FiveLayerReLU(c_in,c_mid,c_mid,c_mid,c_mid,c_mid,1)
    else:
        print('Need to implement. Using a single-layer network')
        net = OneLayerReLU(c_in,c_mid,1)
    optimizer = torch.optim.SGD(net.parameters(), lr=learning_rate, momentum=0.9, weight_decay=weight_decay)
    if opt_adam==True:
        optimizer = torch.optim.Adam(net.parameters(), lr=learning_rate, weight_decay=weight_decay)  
    loss_func = torch.nn.MSELoss()
    for epoch in range(N_epoch):
        perm = torch.randperm(N)
        sum_loss = 0
        for i in range(0, N, batch_size):
            x = X[perm[i : i + batch_size],:]
            y = Y[perm[i : i + batch_size]]
            optimizer.zero_grad()
            output = net(x).squeeze()
            loss = loss_func(output,y)
            loss.backward()
            optimizer.step()
            sum_loss += loss.data.numpy()
        if(print_mid_loss==True):
            print("Epoch:{:4d}\tloss:{}".format(epoch, sum_loss / N))
    cc_sum = 0
    err_sum = 0
    for i in range(0, N, batch_size):
        x = X[perm[i : i + batch_size],:]
        y = Y[perm[i : i + batch_size]]
        output = net(x).squeeze()
        err = (output-y).abs().sum().detach()
        pre_r = torch.round(output)
        cc = (pre_r==y)
        cc_sum += cc.sum().float()
        err_sum += err
    acc=cc_sum/N
    err=err_sum/N
    return net,acc,err


def NN_test(X,Y,net,batch_size=32):
    N = X.shape[0]
    perm = torch.randperm(N)
    cc_sum = 0
    err_sum = 0
    for i in range(0, N, batch_size):
        x = X[perm[i : i + batch_size],:]
        y = Y[perm[i : i + batch_size]]
        output = net(x).squeeze()
        err = (output-y).abs().sum().detach()
        pre_r = torch.round(output)
        cc = (pre_r==y)
        cc_sum += cc.sum().float()
        err_sum += err
    acc=cc_sum/N
    err=err_sum/N
    return acc,err


def SVM_train(X,Y,ls_lambda=0.001,learning_rate=0.001,N_epoch=10,batch_size=128,opt_adam=True,print_mid_loss=False):
    N = X.shape[0]
    c_in = X.shape[1]
    net = LinearF(c_in,1)
    optimizer = torch.optim.SGD(net.parameters(), learning_rate)
    if opt_adam==True:
        optimizer = torch.optim.Adam(net.parameters(), learning_rate)
    for epoch in range(N_epoch):
        perm = torch.randperm(N)
        sum_loss = 0
        for i in range(0, N, batch_size):
            x = X[perm[i : i + batch_size],:]
            y = Y[perm[i : i + batch_size]]
            optimizer.zero_grad()
            #print(x.shape)
            output = net(x)
            #print(output.shape)
            loss = torch.mean(torch.clamp(1 - output.t() * y, min=0))  # hinge loss
            loss += ls_lambda * torch.mean(net.linear1.weight ** 2)  # l2 penalty
            loss.backward()
            optimizer.step()
            sum_loss += loss.data.numpy()
        if(print_mid_loss==True):
            print("Epoch:{:2d}\tloss:{}".format(epoch, sum_loss / N))
    cc_sum = 0
    for i in range(0, N, batch_size):
        x = X[perm[i : i + batch_size],:]
        y = Y[perm[i : i + batch_size]]
        pre = net(x)
        pre_r = (pre>0).squeeze()
        pre_r = ((pre_r.float())-0.5)*2
        cc = (pre_r==y)
        cc_sum += cc.sum().float()
        #print(cc_sum)
    acc = cc_sum.float()/Y.shape[0]
    if(print_mid_loss==True):
        print(net.linear1.weight)
        print(net.linear1.bias)
        print(acc)
    return net,pre,acc

def SVM_train_GPU(X,Y,ls_lambda=0.001,learning_rate=0.001,N_epoch=10,batch_size=128,opt_adam=True,print_mid_loss=False):
    N = X.shape[0]
    c_in = X.shape[1]
    net = LinearF(c_in,1).cuda()
    optimizer = torch.optim.SGD(net.parameters(), learning_rate)
    if opt_adam==True:
        optimizer = torch.optim.Adam(net.parameters(), learning_rate)
    for epoch in range(N_epoch):
        perm = torch.randperm(N)
        sum_loss = 0
        for i in range(0, N, batch_size):
            x = X[perm[i : i + batch_size],:].cuda()
            y = Y[perm[i : i + batch_size]].cuda()
            optimizer.zero_grad()
            #print(x.shape)
            output = net(x)
            #print(output.shape)
            loss = torch.mean(torch.clamp(1 - output.t() * y, min=0))  # hinge loss
            loss += ls_lambda * torch.mean(net.linear1.weight ** 2)  # l2 penalty
            loss.backward()
            optimizer.step()
            sum_loss += loss.data.cpu().numpy()
        if(print_mid_loss==True):
            print("Epoch:{:2d}\tloss:{}".format(epoch, sum_loss / N))
    cc_sum = 0
    for i in range(0, N, batch_size):
        x = X[perm[i : i + batch_size],:].cuda()
        y = Y[perm[i : i + batch_size]].cuda()
        pre = net(x)
        pre_r = (pre>0).squeeze()
        pre_r = ((pre_r.float())-0.5)*2
        cc = (pre_r==y)
        cc_sum += cc.sum().cpu().float()
    acc = cc.sum()/Y.shape[0]
    if(print_mid_loss==True):
        print(net.linear1.weight)
        print(net.linear1.bias)
        print(acc)
    return net,pre,acc


def CLS_N_train_GPU(X,Y,n_class=31,n_layer=1,c_mid=10,learning_rate=0.001,weight_decay=0.01,N_epoch=1,batch_size=128,opt_adam=True,print_mid_loss=False):
    # Binary Classification Neural Network Code
    N = X.shape[0]
    c_in = X.shape[1]
    if n_layer==1:
        net = OneLayerReLU(c_in,c_mid,n_class).cuda()
    elif n_layer==2:
        net = TwoLayerReLU(c_in,c_mid,c_mid,n_class).cuda()
    elif n_layer==3:
        net = ThreeLayerReLU(c_in,c_mid,c_mid,c_mid,n_class).cuda()
    elif n_layer==4:
        net = FourLayerReLU(c_in,c_mid,c_mid,c_mid,c_mid,n_class).cuda()
    elif n_layer==5:
        net = FiveLayerReLU(c_in,c_mid,c_mid,c_mid,c_mid,c_mid,n_class).cuda()
    else:
        print('Need to implement. Using a single-layer network')
        net = OneLayerReLU(c_in,c_mid,n_class).cuda()
    optimizer = torch.optim.SGD(net.parameters(), learning_rate)
    if opt_adam==True:
        optimizer = torch.optim.Adam(net.parameters(), learning_rate)
    loss_function = nn.CrossEntropyLoss()  
    for epoch in range(N_epoch):
        perm = torch.randperm(N)
        sum_loss = 0
        for i in range(0, N, batch_size):
            x = X[perm[i : i + batch_size],:].cuda()
            y = Y[perm[i : i + batch_size]].cuda()
            optimizer.zero_grad()
            #print(x.shape)
            output = net(x)
            #print(output.shape)
            labels = y.reshape([-1]).long()
            loss = loss_function(output, labels)
            loss.backward()
            optimizer.step()
            sum_loss += loss.data.cpu().numpy()
        if(print_mid_loss==True):
            print("Epoch:{:4d}\tloss:{}".format(epoch, sum_loss / N))
    cc_sum = 0
    for i in range(0, N, batch_size):
        x = X[perm[i : i + batch_size],:].cuda()
        y = Y[perm[i : i + batch_size]].cuda()
        pre = net(x)
        labels = y.reshape([-1]).long()
        predicted = torch.argmax(pre.data,1)   #different
        cc = torch.eq(predicted,labels).int().sum().item()
        cc_sum += cc
        #print(cc_sum)
    acc = cc_sum/Y.shape[0]
    if(print_mid_loss==True):
        print(acc)
    return net,acc

def CLS_test_GPU(X,Y,net,batch_size=32):
    N = X.shape[0]
    perm = torch.randperm(N)
    cc_sum = 0
    for i in range(0, N, batch_size):
        x = X[perm[i : i + batch_size],:].cuda()
        y = Y[perm[i : i + batch_size]].cuda()
        pre = net(x)
        labels = y.reshape([-1]).long()
        predicted = torch.argmax(pre.data,1) 
        cc = torch.eq(predicted,labels).int().sum().item()
        cc_sum += cc
    acc=cc_sum/N
    return acc

def CLS_grad_GPU(X,Y,net,batch_size=1024):
    N = X.shape[0]
    perm = torch.randperm(N)
    cc_sum = 0
    optimizer = torch.optim.Adam(net.parameters(), 0.001)
    loss_function = nn.CrossEntropyLoss()
    for i in range(0, N, batch_size):
        x = X[perm[i : i + batch_size],:].cuda()
        y = Y[perm[i : i + batch_size]].cuda()
        x.requires_grad = True
        optimizer.zero_grad()
        output = net(x)
        labels = y.reshape([-1]).long()
        loss = loss_function(output, labels)
        loss.backward()
        if i == 0:
            sen_vec = x.grad.abs().sum(dim=0)
        else:
            sen_vec = sen_vec + x.grad.abs().sum(dim=0)
    sen_vec=sen_vec/N*batch_size
    return sen_vec

def REG_train_GPU(X,Y,n_layer=1,c_mid=10,learning_rate=0.001,weight_decay=0.01,N_epoch=1,batch_size=128,opt_adam=True,print_mid_loss=False):
    # Regression Neural Network Code
    N = X.shape[0]
    c_in = X.shape[1]
    if n_layer==1:
        net = OneLayerReLU(c_in,c_mid,1).cuda()
    elif n_layer==2:
        net = TwoLayerReLU(c_in,c_mid,c_mid,1).cuda()
    elif n_layer==3:
        net = ThreeLayerReLU(c_in,c_mid,c_mid,c_mid,1).cuda()
    elif n_layer==4:
        net = FourLayerReLU(c_in,c_mid,c_mid,c_mid,c_mid,1).cuda()
    elif n_layer==5:
        net = FiveLayerReLU(c_in,c_mid,c_mid,c_mid,c_mid,c_mid,1).cuda()
    else:
        print('Need to implement. Using a single-layer network')
        net = OneLayerReLU(c_in,c_mid,n_class).cuda()
    optimizer = torch.optim.SGD(net.parameters(), learning_rate)
    if opt_adam==True:
        optimizer = torch.optim.Adam(net.parameters(), learning_rate)
    loss_function = torch.nn.MSELoss()
    for epoch in range(N_epoch):
        perm = torch.randperm(N)
        sum_loss = 0
        for i in range(0, N, batch_size):
            x = X[perm[i : i + batch_size],:].cuda()
            y = Y[perm[i : i + batch_size]].cuda()
            optimizer.zero_grad()
            #print(x.shape)
            output = net(x).squeeze()
            #print(output.shape)
            loss = loss_function(output,y)
            loss.backward()
            optimizer.step()
            sum_loss += loss.data.cpu().numpy()
        if(print_mid_loss==True):
            print("Epoch:{:4d}\tloss:{}".format(epoch, sum_loss / N))
    cc_sum = 0
    err_sum = 0
    for i in range(0, N, batch_size):
        x = X[perm[i : i + batch_size],:].cuda()
        y = Y[perm[i : i + batch_size]].cuda()
        output = net(x).squeeze()
        err = (output-y).abs().sum().detach().item()
        pre_r = torch.round(output)
        cc = (pre_r==y)
        cc_sum += cc.sum().float().item()
        err_sum += err
    acc=cc_sum/N
    err=err_sum/N
    if(print_mid_loss==True):
        print(acc)
        print(err)
    return net,acc,err

def REG_test_GPU(X,Y,net,batch_size=32):
    N = X.shape[0]
    perm = torch.randperm(N)
    cc_sum = 0
    err_sum = 0
    for i in range(0, N, batch_size):
        x = X[perm[i : i + batch_size],:].cuda()
        y = Y[perm[i : i + batch_size]].cuda()
        output = net(x).squeeze()
        err = (output-y).abs().sum().detach().item()
        pre_r = torch.round(output)
        cc = (pre_r==y)
        cc_sum += cc.sum().float().item()
        err_sum += err
    acc=cc_sum/N
    err=err_sum/N
    return acc,err        
        

def REG_grad_GPU(X,Y,net,batch_size=1024):
    N = X.shape[0]
    perm = torch.randperm(N)
    optimizer = torch.optim.Adam(net.parameters(), 0.001)
    loss_function = torch.nn.MSELoss()
    for i in range(0, N, batch_size):
        x = X[perm[i : i + batch_size],:].cuda()
        y = Y[perm[i : i + batch_size]].cuda()
        x.requires_grad = True
        optimizer.zero_grad()
        output = net(x)
        obj_f = output.sum()
        obj_f.backward()
        if i == 0:
            sen_vec = x.grad.abs().sum(dim=0)
        else:
            sen_vec = sen_vec + x.grad.abs().sum(dim=0)
    sen_vec=sen_vec/N
    return sen_vec

if __name__ == "__main__":
    # Test the square_root function
    print('Functions for the Mathine Learning module')
