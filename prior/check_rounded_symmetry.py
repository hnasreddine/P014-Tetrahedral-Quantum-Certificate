from pathlib import Path
p=Path(__file__).parent/'rounded_grams.txt'
s=p.read_text().split();j=1
for _ in range(int(s[0])):
 name=s[j];n=int(s[j+1]);j+=2
 m=list(map(int,s[j:j+n*n]));j+=n*n;r=n//2
 assert all(m[i*n+k]==m[k*n+i] for i in range(n) for k in range(n))
 assert all(m[i*n+k]==m[(i+r)*n+k+r] and m[i*n+k+r]==-m[(i+r)*n+k] for i in range(r) for k in range(r))
 print(name,'EXACT REAL-HERMITIAN-EMBEDDING PASS')
