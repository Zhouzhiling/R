���ڴ���Ĵ�����Cluster.r��������Ҫ������fun.Calculate <- function(TextNum,KeyWord,CellThresholdUpper,CellThresholdLower,GeneThreshold,ClusterNum)

#����һn�����ݵ�ǰn�����ı���Ϣ������ʱ��Ҫɾ��
#������m����������ÿ���������������Ҫ��Ψһ�ģ�
#������k1��ȥ��������ϸ�� total reads < k
#������k2��ȥ��������ϸ�� total reads > k
#������p��ȥ���ͱ����� total count < p
#������q��ѡ�񷽲�Ϊǰq�Ļ�����о���


100��101��111��129��157�е�·����Ҫ������ǰ�޸�

����Ľ��������ļ�ClusterResult.csv��

DESeq�Ĵ��������ļ�DESeqResult.csv��

�����ڵ������ǣ�
��֪��SC3����Ľ����ֱ�����������и�ʽҪ��Ļ���Ҫ��111�������
