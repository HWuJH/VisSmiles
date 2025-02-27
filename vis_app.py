import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
import py3Dmol  # 兼容 Streamlit Cloud

# 设置页面标题
st.title("🔬 SMILES 处理工具")

# 侧边栏输入
st.sidebar.header("📝 输入 SMILES")
smiles = st.sidebar.text_input("请输入 SMILES:", value="CCO")  # 默认乙醇 (CCO)

# 获取分子对象的函数
def get_mol(smiles):
    """将 SMILES 转换为 RDKit 分子对象"""
    try:
        return Chem.MolFromSmiles(smiles) if smiles else None
    except:
        return None

# 处理 SMILES
mol = get_mol(smiles)

# 侧边栏操作区
st.sidebar.header("📌 操作")

if st.sidebar.button("规范化 SMILES"):
    if mol:
        st.session_state["canonical_smiles"] = Chem.MolToSmiles(mol, canonical=True)
    else:
        st.session_state["canonical_smiles"] = "❌ 无效的 SMILES"

if st.sidebar.button("显示 2D 结构"):
    if mol:
        img = Draw.MolToImage(mol, size=(300, 300))
        st.session_state["mol_2d"] = img
    else:
        st.session_state["mol_2d"] = None

if st.sidebar.button("显示 3D 结构"):
    if mol:
        try:
            mol_3d = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol_3d, AllChem.ETKDG())  # 生成 3D 坐标
            mol_block = Chem.MolToMolBlock(mol_3d)
            if mol_block:  # 确保 mol_block 不为空
                # 生成 HTML + JS 代码
                html_3d = f"""
                <div style="width: 100%; height: 400px;" id="mol3d"></div>
                <script src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>
                <script>
                    var viewer = $3Dmol.createViewer("mol3d", {{ backgroundColor: "white" }});
                    viewer.addModel(`{mol_block}`, "mol");
                    viewer.setStyle({{stick: {}}});
                    viewer.zoomTo();
                    viewer.render();
                </script>
                """
                st.session_state["mol_3d"] = html_3d
            else:
                st.session_state["mol_3d"] = "⚠️ 生成 3D 坐标失败"
        except Exception as e:
            st.session_state["mol_3d"] = f"⚠️ 3D 可视化失败: {e}"

# **调整分区布局**
col1, col2, col3 = st.columns([1.2, 1, 1.5])  # 让 3D 结构区域更大

# **规范化 SMILES 显示**
with col1:
    st.subheader("✅ 规范化 SMILES")
    if "canonical_smiles" in st.session_state:
        st.code(st.session_state["canonical_smiles"], language="markdown")

# **2D 结构显示**
with col2:
    st.subheader("🧪 2D 结构")
    if "mol_2d" in st.session_state and st.session_state["mol_2d"]:
        st.image(st.session_state["mol_2d"], caption="2D 结构", use_column_width=True)

# **3D 结构显示**
with col3:
    st.subheader("🧩 3D 结构")
    if "mol_3d" in st.session_state:
        if "⚠️" in st.session_state["mol_3d"]:  # 处理错误信息
            st.error(st.session_state["mol_3d"])
        else:
            st.components.v1.html(st.session_state["mol_3d"], height=400, scrolling=False)

