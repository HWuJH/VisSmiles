import streamlit as st
from rdkit import Chem
from rdkit.Chem import AllChem
import py3Dmol
from stmol import showmol
import pandas as pd
from io import StringIO

# 设置页面标题
st.set_page_config(page_title="SMILES 可视化与规范化工具", layout="wide")
st.title("🔬 SMILES 可视化与规范化工具")

# 介绍
st.markdown("""
**功能说明：**
1️⃣ **可视化**：将 SMILES 解析为分子图像  
2️⃣ **规范化**：转换为标准化的 SMILES  
3️⃣ **3D 结构生成**：生成 3D 分子模型，并可优化构象  
4️⃣ **文件上传**：批量处理 SMILES 并下载规范化结果  
""")

# 初始化 session_state
if "normalized_smiles" not in st.session_state:
    st.session_state["normalized_smiles"] = []

# **输入框**
smiles = st.text_input("🧪 请输入 SMILES 字符串:", value="CCO")

# **创建功能分区**
col1, col2, col3 = st.columns(3)

# **可视化 SMILES**
with col1:
    if st.button("👁 可视化分子"):
        if not smiles.strip():
            st.error("❌ 请输入有效的 SMILES！")
        else:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                # 生成 3D 结构（伪 2D 视图）
                AllChem.EmbedMolecule(mol)
                mol_block = Chem.MolToMolBlock(mol)
                viewer = py3Dmol.view(width=400, height=400)
                viewer.addModel(mol_block, "mol")
                viewer.setStyle({"stick": {}})
                viewer.setViewStyle({"style": "stick", "rotation": 90})  # 旋转以接近 2D 视图
                viewer.zoomTo()
                showmol(viewer, width=400, height=400)
            else:
                st.error("❌ 无效的 SMILES，请检查格式！")

# **规范化 SMILES**
with col2:
    if st.button("🔄 规范化 SMILES"):
        if not smiles.strip():
            st.error("❌ 请输入有效的 SMILES！")
        else:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                canonical_smiles = Chem.MolToSmiles(mol, canonical=True)
                st.session_state["normalized_smiles"].append([smiles, canonical_smiles])
                st.success(f"✅ 规范化后的 SMILES: `{canonical_smiles}`")
            else:
                st.error("❌ 无效的 SMILES，请重新输入！")

# **生成 3D 结构**
with col3:
    if st.button("🧩 生成 3D 结构"):
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            try:
                AllChem.EmbedMolecule(mol)
                AllChem.UFFOptimizeMolecule(mol)  # 进行 3D 优化
                mol_block = Chem.MolToMolBlock(mol)
                viewer = py3Dmol.view(width=400, height=400)
                viewer.addModel(mol_block, "mol")
                viewer.setStyle({"stick": {}})
                viewer.zoomTo()
                showmol(viewer, width=400, height=400)
            except Exception as e:
                st.error(f"❌ 生成失败：{e}")
        else:
            st.error("❌ 无效的 SMILES！")

# **文件上传功能**
st.sidebar.header("📂 批量处理 SMILES")
uploaded_file = st.sidebar.file_uploader("📥 上传 CSV/TXT 文件 (每行一个 SMILES)", type=["csv", "txt"])

if uploaded_file:
    try:
        smiles_list = uploaded_file.read().decode("utf-8").splitlines()
        st.sidebar.success(f"✅ 读取 {len(smiles_list)} 条 SMILES 数据")
        normalized_data = []

        for smi in smiles_list:
            mol = Chem.MolFromSmiles(smi)
            if mol:
                normalized_data.append([smi, Chem.MolToSmiles(mol, canonical=True)])
            else:
                normalized_data.append([smi, "❌ 无效 SMILES"])

        # 保存到 session_state 以便下载
        st.session_state["normalized_smiles"] = normalized_data

        # 显示数据
        df = pd.DataFrame(normalized_data, columns=["原始 SMILES", "规范化 SMILES"])
        st.sidebar.write(df)

    except Exception as e:
        st.sidebar.error(f"❌ 处理失败：{e}")

# **下载按钮**
if st.sidebar.button("📥 下载规范化结果"):
    if st.session_state["normalized_smiles"]:
        df = pd.DataFrame(st.session_state["normalized_smiles"], columns=["原始 SMILES", "规范化 SMILES"])
        csv = df.to_csv(index=False)
        st.sidebar.download_button("⬇️ 点击下载", csv, "normalized_smiles.csv", "text/csv")
    else:
        st.sidebar.warning("❗ 没有数据可下载！")

# **版权信息**
st.markdown("---")
st.markdown("📌 **开发者：你的名字**  |  💡 **使用 RDKit & Py3DMol**")
