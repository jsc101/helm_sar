' helm_builder.bas
' Import into the Excel workbook via VBE (Alt+F11 -> File -> Import File)
' Assign BuildHELM to a button on the Design sheet.
Option Explicit

' ── Public entry point ───────────────────────────────────────────────────────

Public Sub BuildHELM()
    Dim wsDesign As Worksheet, wsHelm As Worksheet
    Dim polyTypes As Object
    Set polyTypes = LoadMonomerDB()

    Set wsDesign = ThisWorkbook.Sheets("Design")
    Set wsHelm   = ThisWorkbook.Sheets("HELM")

    ' Clear HELM sheet below header
    If wsHelm.Cells(wsHelm.Rows.Count, 1).End(xlUp).Row > 1 Then
        wsHelm.Rows("2:" & wsHelm.Rows.Count).ClearContents
    End If

    ' Read Design sheet headers from row 1
    Dim lastCol As Long
    lastCol = wsDesign.Cells(1, wsDesign.Columns.Count).End(xlToLeft).Column
    Dim headers() As String
    ReDim headers(1 To lastCol) As String
    Dim c As Long
    For c = 1 To lastCol
        headers(c) = Trim(CStr(wsDesign.Cells(1, c).Value))
    Next c

    ' Process each data row
    Dim r As Long, outRow As Long
    outRow = 2
    r = 2
    Do While Trim(CStr(wsDesign.Cells(r, 1).Value)) <> ""
        Dim rowName As String
        rowName = CStr(wsDesign.Cells(r, 1).Value)

        Dim rowData As Object
        Set rowData = CreateObject("Scripting.Dictionary")
        For c = 1 To lastCol
            rowData(headers(c)) = Trim(CStr(wsDesign.Cells(r, c).Value))
        Next c

        Dim helm As String, status As String
        helm   = GenerateHELM(rowData, polyTypes)
        status = ValidateHELM(helm)

        wsHelm.Cells(outRow, 1).Value = rowName
        wsHelm.Cells(outRow, 2).Value = helm
        wsHelm.Cells(outRow, 3).Value = status
        outRow = outRow + 1
        r = r + 1
    Loop

    MsgBox "Done. " & (outRow - 2) & " compounds written to HELM sheet.", vbInformation
End Sub

' ── MonomerDB ────────────────────────────────────────────────────────────────

Private Function LoadMonomerDB() As Object
    Dim ws As Worksheet
    Dim db As Object
    Set db = CreateObject("Scripting.Dictionary")

    On Error Resume Next
    Set ws = ThisWorkbook.Sheets("MonomerDB")
    On Error GoTo 0
    If ws Is Nothing Then
        MsgBox "MonomerDB sheet not found. Run make_template.py to regenerate.", vbExclamation
        Set LoadMonomerDB = db
        Exit Function
    End If

    Dim lastR As Long
    lastR = ws.Cells(ws.Rows.Count, 1).End(xlUp).Row
    Dim i As Long
    For i = 2 To lastR
        Dim sym As String, pt As String
        sym = Trim(CStr(ws.Cells(i, 1).Value))
        pt  = Trim(CStr(ws.Cells(i, 2).Value))
        If sym <> "" Then db(sym) = pt
    Next i
    Set LoadMonomerDB = db
End Function

' ── Row -> HELM ──────────────────────────────────────────────────────────────

Private Function GenerateHELM(rowData As Object, polyTypes As Object) As String
    ' 1. Build main chain
    Dim mainChain As String
    mainChain = ""
    Dim pos As Long
    For pos = 1 To 200
        Dim posKey As String
        posKey = CStr(pos)
        If Not rowData.Exists(posKey) Then Exit For
        Dim aa As String
        aa = rowData(posKey)
        If aa = "" Then Exit For
        If Len(mainChain) > 0 Then mainChain = mainChain & "."
        mainChain = mainChain & FmtSym(aa)
    Next pos

    Dim chains As String, conns As String
    chains = "PEPTIDE1{" & mainChain & "}"
    conns  = ""

    Dim peptideN As Long, chemN As Long
    peptideN = 2
    chemN = 1

    ' 2. Process up to 3 sidechain blocks
    Dim scIdx As Long
    For scIdx = 1 To 3
        Dim sfx As String
        If scIdx = 1 Then
            sfx = ""
        Else
            sfx = "_" & scIdx
        End If

        Dim siteKey As String
        siteKey = "Site" & sfx
        If Not rowData.Exists(siteKey) Then GoTo NextSC
        Dim siteVal As String
        siteVal = rowData(siteKey)
        If siteVal = "" Then GoTo NextSC
        Dim site As Long
        site = CLng(siteVal)

        ' Collect bonds and monomers
        Dim bonds(1 To 20) As String
        Dim monomers(1 To 20) As String
        Dim nSC As Long
        nSC = 0
        Dim k As Long
        For k = 1 To 20
            Dim bKey As String, mKey As String
            bKey = "b" & k & sfx
            mKey = "SC" & k & sfx
            If Not rowData.Exists(bKey) Then Exit For
            Dim bVal As String, mVal As String
            bVal = rowData(bKey)
            mVal = rowData(mKey)
            If bVal = "" Or mVal = "" Then Exit For
            nSC = nSC + 1
            bonds(nSC) = bVal
            monomers(nSC) = mVal
        Next k

        If nSC = 0 Then GoTo NextSC

        If CanCollapse(bonds, monomers, nSC, polyTypes) Then
            Dim scBody As String
            scBody = ""
            For k = 1 To nSC
                If Len(scBody) > 0 Then scBody = scBody & "."
                scBody = scBody & FmtSym(monomers(k))
            Next k
            Dim chainId As String
            chainId = "PEPTIDE" & peptideN
            peptideN = peptideN + 1
            chains = chains & "|" & chainId & "{" & scBody & "}"
            Dim b1Parts() As String
            b1Parts = Split(bonds(1), "-")
            Dim proxRg As String
            proxRg = b1Parts(0)
            If Len(conns) > 0 Then conns = conns & "|"
            conns = conns & chainId & ",PEPTIDE1,1:R1-" & site & ":" & proxRg
        Else
            ' Expand: one chain per monomer
            Dim prevId As String
            prevId = ""
            Dim isFirst As Boolean
            isFirst = True

            For k = 1 To nSC
                Dim ptype As String
                ptype = "PEPTIDE"
                If polyTypes.Exists(monomers(k)) Then ptype = polyTypes(monomers(k))

                Dim curId As String
                If ptype = "CHEM" Then
                    curId = "CHEM" & chemN
                    chemN = chemN + 1
                Else
                    curId = "PEPTIDE" & peptideN
                    peptideN = peptideN + 1
                End If

                chains = chains & "|" & curId & "{" & FmtSym(monomers(k)) & "}"

                Dim bParts() As String
                bParts = Split(bonds(k), "-")
                Dim pRg As String, dRg As String
                pRg = bParts(0)
                dRg = bParts(1)

                Dim conn As String
                If isFirst Then
                    conn = curId & ",PEPTIDE1,1:" & dRg & "-" & site & ":" & pRg
                    isFirst = False
                Else
                    conn = curId & "," & prevId & ",1:" & dRg & "-1:" & pRg
                End If
                If Len(conns) > 0 Then conns = conns & "|"
                conns = conns & conn
                prevId = curId
            Next k
        End If

NextSC:
    Next scIdx

    If Len(conns) > 0 Then
        GenerateHELM = chains & "$" & conns & "$$$V2.0"
    Else
        GenerateHELM = chains & "$$$$V2.0"
    End If
End Function

' ── Helpers ──────────────────────────────────────────────────────────────────

Private Function FmtSym(sym As String) As String
    If Len(sym) = 1 Then
        FmtSym = sym
    Else
        FmtSym = "[" & sym & "]"
    End If
End Function

Private Function CanCollapse(bonds() As String, monomers() As String, _
                              nSC As Long, polyTypes As Object) As Boolean
    Dim k As Long
    For k = 1 To nSC
        Dim pt As String
        pt = "PEPTIDE"
        If polyTypes.Exists(monomers(k)) Then pt = polyTypes(monomers(k))
        If pt <> "PEPTIDE" Then
            CanCollapse = False
            Exit Function
        End If
    Next k
    ' b1 distal Rg must be R1
    Dim b1p() As String
    b1p = Split(bonds(1), "-")
    If b1p(1) <> "R1" Then
        CanCollapse = False
        Exit Function
    End If
    ' All subsequent bonds must be R2-R1
    For k = 2 To nSC
        If bonds(k) <> "R2-R1" Then
            CanCollapse = False
            Exit Function
        End If
    Next k
    CanCollapse = True
End Function

Private Function ValidateHELM(helm As String) As String
    ' Lightweight check: V2.0 suffix and balanced braces
    If Right(helm, 6) <> "V2.0" Then
        ValidateHELM = "ERROR: missing V2.0 suffix"
        Exit Function
    End If
    Dim depth As Long
    depth = 0
    Dim i As Long
    For i = 1 To Len(helm)
        Dim ch As String
        ch = Mid(helm, i, 1)
        If ch = "{" Then depth = depth + 1
        If ch = "}" Then depth = depth - 1
        If depth < 0 Then
            ValidateHELM = "ERROR: unbalanced braces"
            Exit Function
        End If
    Next i
    If depth <> 0 Then
        ValidateHELM = "ERROR: unbalanced braces"
    Else
        ValidateHELM = "OK"
    End If
End Function
