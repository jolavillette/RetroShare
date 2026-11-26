/*******************************************************************************
 * retroshare-gui/src/gui/chat/ChatLobbyPolicyDialog.cpp                       *
 *                                                                             *
 * Copyright (C) 2026, Retroshare Team <retroshare.project@gmail.com>          *
 *                                                                             *
 * This program is free software: you can redistribute it and/or modify        *
 * it under the terms of the GNU Affero General Public License as              *
 * published by the Free Software Foundation, either version 3 of the          *
 * License, or (at your option) any later version.                             *
 *                                                                             *
 * This program is distributed in the hope that it will be useful,             *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of              *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the                *
 * GNU Affero General Public License for more details.                         *
 *                                                                             *
 * You should have received a copy of the GNU Affero General Public License    *
 * along with this program. If not, see <https://www.gnu.org/licenses/>.       *
 *                                                                             *
 *******************************************************************************/
#include <QCheckBox>
#include <QDialogButtonBox>
#include <QHeaderView>
#include <QLabel>
#include <QTreeWidget>
#include <QVBoxLayout>

#include <retroshare/rsidentity.h>
#include <retroshare/rspeers.h>

#include "ChatLobbyPolicyDialog.h"

#define COLUMN_NODE     0
#define COLUMN_IDENTITY 1
#define COLUMN_ADMIN    2

#define ROLE_NODE_ID    Qt::UserRole

ChatLobbyPolicyDialog::ChatLobbyPolicyDialog(const ChatLobbyId& lobbyId, QWidget *parent)
	: QDialog(parent, Qt::WindowSystemMenuHint | Qt::WindowTitleHint | Qt::WindowCloseButtonHint),
	  mLobbyId(lobbyId), mEnabled(nullptr), mMembersSee(nullptr), mMembersTalk(nullptr), mParticipants(nullptr)
{
	setWindowTitle(tr("Hosted room"));
	new QVBoxLayout(this);

	ChatLobbyRole role;

	// A node that holds rules of its own edits them, whatever it was told.
	if(!rsChats->getChatLobbyPolicy(mLobbyId, mPolicy) && rsChats->getChatLobbyRole(mLobbyId, role) && role.ruled)
		buildRoleView(role);
	else
		buildEditor();
}

void ChatLobbyPolicyDialog::buildRoleView(const ChatLobbyRole& role)
{
	QStringList lines;
	lines << tr("The node %1 is the host of this room: it decides what it passes on between the participants.").arg(QString::fromUtf8(rsPeers->getPeerName(role.host).c_str()));

	if(role.admin)
		lines << tr("You are an admin: what you write goes to everyone, and you receive everything.");
	else
	{
		lines << tr("You are a member.");
		lines << (role.membersTalk ? tr("What you write goes to everyone.") : tr("What you write reaches the admins only."));
		lines << (role.membersSeeEachOther ? tr("You see every participant.") : tr("You see the admins only."));
	}

	QLabel *label = new QLabel(lines.join("\n\n"));
	label->setWordWrap(true);
	layout()->addWidget(label);

	QDialogButtonBox *buttons = new QDialogButtonBox(QDialogButtonBox::Close);
	connect(buttons, SIGNAL(rejected()), this, SLOT(reject()));
	layout()->addWidget(buttons);
}

void ChatLobbyPolicyDialog::buildEditor()
{
	QLabel *label = new QLabel(tr("As the host of this room, your node decides what it passes on between the participants. What an admin writes goes to everyone, and admins receive everything; the host is always an admin. The other participants are members.\n\nThis holds for participants that reach the room through your node. Participants connected to each other directly exchange everything."));
	label->setWordWrap(true);
	layout()->addWidget(label);

	mEnabled = new QCheckBox(tr("My node is the host of this room"));
	mMembersSee = new QCheckBox(tr("Members see each other"));
	mMembersTalk = new QCheckBox(tr("Members talk to each other"));

	mMembersSee->setToolTip(tr("If unchecked, a member only sees the admins in the participant list."));
	mMembersTalk->setToolTip(tr("If unchecked, what a member writes reaches the admins only."));

	mEnabled->setChecked(mPolicy.enabled);
	mMembersSee->setChecked(mPolicy.membersSeeEachOther);
	mMembersTalk->setChecked(mPolicy.membersTalk);

	layout()->addWidget(mEnabled);
	layout()->addWidget(mMembersSee);
	layout()->addWidget(mMembersTalk);

	mParticipants = new QTreeWidget;
	mParticipants->setRootIsDecorated(false);
	mParticipants->setHeaderLabels(QStringList() << tr("Invited node") << tr("Identity in the room") << tr("Admin"));
	mParticipants->header()->setSectionResizeMode(COLUMN_NODE, QHeaderView::Stretch);
	mParticipants->header()->setSectionResizeMode(COLUMN_IDENTITY, QHeaderView::Stretch);
	mParticipants->header()->setStretchLastSection(false);
	layout()->addWidget(mParticipants);

	fillParticipants(mPolicy);

	QDialogButtonBox *buttons = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel);
	connect(buttons, SIGNAL(accepted()), this, SLOT(apply()));
	connect(buttons, SIGNAL(rejected()), this, SLOT(reject()));
	layout()->addWidget(buttons);

	connect(mEnabled, SIGNAL(toggled(bool)), this, SLOT(enableChanged(bool)));
	enableChanged(mEnabled->isChecked());

	resize(560, 420);
}

void ChatLobbyPolicyDialog::fillParticipants(const ChatLobbyPolicy& policy)
{
	for(const ChatLobbyParticipant& participant : policy.participants)
	{
		QTreeWidgetItem *item = new QTreeWidgetItem(mParticipants);

		RsPeerDetails node;
		QString nodeName = QString::fromStdString(participant.node.toStdString());

		if(rsPeers->getPeerDetails(participant.node, node))
			nodeName = QString::fromUtf8(node.name.c_str()) + " (" + QString::fromUtf8(node.location.c_str()) + ")";

		QString identity = tr("not known yet");
		RsIdentityDetails details;

		if(!participant.identity.isNull())
			identity = rsIdentity->getIdDetails(participant.identity, details)
			         ? QString::fromUtf8(details.mNickname.c_str())
			         : QString::fromStdString(participant.identity.toStdString());

		item->setText(COLUMN_NODE, nodeName);
		item->setData(COLUMN_NODE, ROLE_NODE_ID, QString::fromStdString(participant.node.toStdString()));
		item->setText(COLUMN_IDENTITY, identity);
		item->setToolTip(COLUMN_IDENTITY, participant.identity.isNull()
		                 ? tr("What this participant writes is passed on as written by a member until its identity in the room is known.")
		                 : QString::fromStdString(participant.identity.toStdString()));
		item->setFlags(item->flags() | Qt::ItemIsUserCheckable);
		item->setCheckState(COLUMN_ADMIN, participant.admin ? Qt::Checked : Qt::Unchecked);
	}
}

void ChatLobbyPolicyDialog::enableChanged(bool enabled)
{
	mMembersSee->setEnabled(enabled);
	mMembersTalk->setEnabled(enabled);
	mParticipants->setEnabled(enabled);
}

void ChatLobbyPolicyDialog::apply()
{
	const bool held = rsChats->getChatLobbyPolicy(mLobbyId, mPolicy);

	// Nothing to hold for a room that never had rules and gets none.
	if(held || mEnabled->isChecked())
	{
		rsChats->setChatLobbyPolicy(mLobbyId, mEnabled->isChecked(), mMembersSee->isChecked(), mMembersTalk->isChecked());
		rsChats->getChatLobbyPolicy(mLobbyId, mPolicy);

		for(int i = 0; i < mParticipants->topLevelItemCount(); ++i)
		{
			QTreeWidgetItem *item = mParticipants->topLevelItem(i);

			const RsPeerId node(item->data(COLUMN_NODE, ROLE_NODE_ID).toString().toStdString());
			const bool admin = item->checkState(COLUMN_ADMIN) == Qt::Checked;

			for(const ChatLobbyParticipant& participant : mPolicy.participants)
				if(participant.node == node && participant.admin != admin)
					rsChats->setChatLobbyParticipant(mLobbyId, node, RsGxsId(), admin);
		}
	}
	accept();
}
