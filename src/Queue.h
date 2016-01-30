# if !defined ___template___CQueue_h___
# define ___template___CQueue_h___

#include <stdio.h>

#define TRUE 1
#define FALSE 0
#define OK 1
#define ERROR 0
#define INFEASIBLE -1
#define OVER_FLOW -2
typedef int Status;


template<class QType>
class CQueue{
public:
        CQueue(unsigned size);
	CQueue();
	~CQueue();
	void EnQueue(QType e); //element will enter queue even if it is full
	void InitQueue(unsigned size);
	void DeleQueue();
	void DeQueue_();
	QType DeQueue();       //element will enter queue even if it is empty
	bool IsEmpty();
	bool IsFull();
	void ViewQueue();
private:
	QType *queue; 
	unsigned front;
	unsigned rear;
	unsigned num; 
};


template<class QType>
CQueue<QType>::CQueue(unsigned size){
  queue = new QType [size];
  front = 0; 
  rear  = 0; 
  num   = size;
}

template<class QType>
CQueue<QType>::CQueue(){
  queue = NULL;
  front = 0; 
  rear  = 0; 
  num   = 0;
}

template<class QType>
void CQueue<QType>::InitQueue(unsigned size){
  queue = new QType [size];
  front = 0;
  rear  = 0;
  num   = size;
}

template<class QType>
void CQueue<QType>::DeleQueue(){
  delete[] queue;
}

template<class QType>
CQueue<QType>::~CQueue(){
}




template <class QType>
void CQueue<QType>::DeQueue_(){
  if(front==rear) return;
  front = (front + 1)%num; 
}

template <class QType>
void CQueue<QType>::EnQueue(QType e){
  queue[rear] = e; 
  rear = (rear + 1)%num; 
}

template <class QType>
QType CQueue<QType>::DeQueue(){
  QType e = queue[front];
  front = (front + 1)%num; 
  return e; 
}

template <class QType>
bool CQueue<QType>::IsEmpty(){
	return front==rear;
}

template <class QType>
bool CQueue<QType>::IsFull(){
	return front==rear;
}

template<class QType>
void CQueue<QType>::ViewQueue(){
  if(0==num)
    return;
  int i=front;
  while(i!=rear){
    printf("%4d ",queue[i]);
    i=(i+1)%num;
  }
  printf("\n");
}

#endif
