##
linux下python调用C接口，SWIG版本
##
步骤：
1.linux下载swig命令  apt install swig
2.新建*.i文件 文件示例见align.i
3.构建so动态库 依次运行命令（1）swig -c++/c -python *.i （2）g++ -fPIC -shared align_wrap.cxx/c(swig生成的和python交互的接口) align.cpp/c(自己写的供python调用的c接口) -o _align.so(你的.so文件名) -I /usr/include/python版本(安装的python目录)
4.将so动态库文件导入系统环境变量 在~/.bashrc文件添加 export LD_LABRARY_PATH=$LD_LABRARY_PATH:你的so文件目录
5.source ~/.bashrc 更新系统环境变量
6.然后python即可调用相应模块模块



