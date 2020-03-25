import numpy as np
import matplotlib.pyplot as plt

def read_fasta(filename):
    s = ''
    for l in open(filename, 'r'):
        if l[0] != '>': s += l.rstrip('\n')
    return s


def dotmat(s,t,W=1,S=1):
    """ Dot Matrix of sequence comparison
    s,t = biological sequence (protein or DNA)
    W = window size; S = stringency
    reference: D. W. Mount "Bioinformatics" p59
    """
    s = np.array(list(s)).reshape(1,-1)
    t = np.array(list(t)).reshape(1,-1)
    for _ in range(W):
        s = np.vstack((s, np.roll(s[-1],-1)))
        t = np.vstack((t, np.roll(t[-1],-1)))

    s = s.T[:-W]
    t = t.T[:-W]

    i,j = np.meshgrid(range(len(s)), range(len(t)), indexing='ij')
    return np.count_nonzero(s[i]==t[j], axis=2) >= S


s = read_fasta('CAA35120.1.fasta')
t = read_fasta('P17678.1.fasta')

plt.spy(dotmat(t,s,3,2))
plt.title('erythroid transcription factor', fontsize=20)
plt.xlabel('homo sapiens', fontsize=20)
plt.ylabel('chicken', fontsize=20)
plt.tight_layout()
plt.savefig('fig1.eps')
plt.show()
