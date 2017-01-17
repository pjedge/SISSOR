coverage_cut = 2
base_letters = 'ACGTacgt'
numbers = '0123456789'
bases = ['A','T','G','C']

import re

caret_money = re.compile(r'\^.|\$')
comma_or_period = re.compile(r'\.|\,')
plus_or_minus = re.compile(r'\+|-')

def parse_mpileup_base_qual(raw_bd, raw_qd, ref_base):

    bd = [] # base data
    qd = [] # qual data

    indel_count = len(re.findall(caret_money, raw_bd))

    if not re.search(plus_or_minus,raw_bd):

        bd = str.upper(re.sub(caret_money,'', raw_bd))
        bd = re.sub(comma_or_period, ref_base, bd)
        qd = [10**((ord(q) - 33) * -0.1) for q in raw_qd]
        paired_bd_qd = [(b,q) for b,q in zip(bd,qd) if b not in ['>','<','*','n','N']]
        if len(paired_bd_qd) < coverage_cut:
            return [],[], 0

        bd, qd = zip(*paired_bd_qd)
        bd = list(bd)
        qd = list(qd)

    else:

        i = 0   # index for base data
        j = 0   # index for qual data

        while i < len(raw_bd):
            if raw_bd[i] == '.' or raw_bd[i] == ',':   # reference base call
                bd.append(ref_base)
                qd.append(10**((ord(raw_qd[j]) - 33) * -0.1))
                i += 1
                j += 1
            elif raw_bd[i] in base_letters:            # variant call
                bd.append(str.upper(raw_bd[i]))
                qd.append(10**((ord(raw_qd[j]) - 33) * -0.1))
                i += 1
                j += 1
            elif raw_bd[i] == '+' or raw_bd[i] == '-': # indel
                indel_count += 1
                num = int(raw_bd[i+1])
                i += 2
                while(raw_bd[i] in numbers):
                    num *= 10
                    num += int(raw_bd[i])
                    i += 1
                i += num
            elif raw_bd[i] == '^':
                i += 2
            elif raw_bd[i] == '$':
                i += 1
            elif raw_bd[i] in '><*nN':                 # reference skip or deletion
                i += 1
                j += 1

        assert(i == len(raw_bd))
        assert(j == len(raw_qd))

    assert(len(bd) == len(qd))

    for b in bd:
        assert(b in bases)

    return bd, qd, indel_count
