# Criar e mudar para o branch dev
git checkout -b dev

# enviar o novo branch para github
git push -u origin dev

# ver os branchs
git branch

# Alternar entre os branchs
git checkout main
# ou
git checkout dev

# COmando mais moderno
git switch main
# ou
git switch dev

####### FAZER PUSH NO GIT HUB
# PUSH NO DEV (so se for para desenvolvimento)
# Correcoes do worflow padrao jogamos no main mesmo

# 1. Faça commit e push no main

git checkout main
git add .
git commit -m "update"
git push origin main

# 2. Vá para o branch dev e atualize-o com o main
# Ou o contrario caso seja novas atualizacoes

git checkout dev
git merge main --no-edit  # Mescla as mudanças do main no dev
git push origin dev  # Envia as mudanças para o remoto

#### PULL FROM ALL BRANCHS
git fetch --all --prune # atualiza os branchs online com o local
git pull --all


